#' Technical and efficiency change estimation
#' 
#' Estimates technical and efficiency change using SNFA
#' 
#' @param df Data frame with variables used in estimation
#' @param input.var.names Names of input variables; must appear in df
#' @param output.var.name Name of output variable; must appear in df
#' @param firm.var.name Name of firm variable; must appear in df
#' @param time.var.name Name of time variable; must appear in df
#' @param method Constraints to apply; "u" for unconstrained, "m" for monotonically increasing, and "mc" for monotonically increasing and concave
#' 
#' @return Returns a data.frame with the following columns
#' \item{firm.var.name}{Column of firm name data}
#' \item{time.var.name}{Column of time period data}
#' \item{efficiency.change}{Average annual efficiency change since the previous period in data}
#' \item{technical.change}{Average annual technical change since the previous period in data}
#' \item{productivity.change}{Average annual productivity change since the previous period in data}
#' 
#' @details 
#' This function decomposes change in productivity into efficiency and technical change, as in
#' Fare et al. (1994), using smooth non-parametric frontier analysis. Denoting \eqn{D_s(x_t, y_t)}
#' as the efficiency of the production plan in year t relative to the production frontier in year s,
#' efficiency change for a given firm in year t is calculated as
#' 
#' \deqn{\frac{D_{t+1}(x_{t+1}, y_{t+1})}{D_t(x_t, y_t)},}
#' 
#' and technical change is given by
#' 
#' \deqn{\left(
#' \frac{D_t(x_{t+1}, y_{t+1})}{D_{t+1}(x_{t+1}, y_{t+1})}\times
#' \frac{D_t(x_t, y_t)}{D_{t+1}(x_t, y_t)}
#' \right)^{1/2}.}
#' 
#' @examples
#' data(panel.production)
#' 
#' results.df <- technical.efficiency.change(df = panel.production,
#'                                           input.var.names = c("X.1", "X.2", "X.3"),
#'                                           output.var.name = "y",
#'                                           firm.var.name = "Firm",
#'                                           time.var.name = "Year")
#' 
#' #Plot changes over time by firm
#' library(ggplot2)
#' 
#' ggplot(results.df, aes(Year, technical.change)) +
#'   geom_line(aes(color = Firm))
#' ggplot(results.df, aes(Year, efficiency.change)) +
#'   geom_line(aes(color = Firm))
#' ggplot(results.df, aes(Year, productivity.change)) +
#'   geom_line(aes(color = Firm))
#'   
#' @references
#' \insertRef{Fare}{snfa}
#' 
#' @export

technical.efficiency.change <-
  function(df, input.var.names, output.var.name, firm.var.name, time.var.name, method = "u"){
    if (class(df) != "data.frame") stop("df must be of class data.frame.")
    if (!(all(input.var.names %in% names(df)))) stop("All input.var.names must be in df.")
    if (!(output.var.name %in% names(df))) stop("output.var.name must be in df.")
    if (!(firm.var.name %in% names(df))) stop("firm.var.name must be in df.")
    if (!(time.var.name %in% names(df))) stop("time.var.name must be in df.")
    
    years <- sort(unique(df[,time.var.name]))
    num.years <- length(years)
    if (num.years <= 1) stop("Must be at least 2 years in df.")
    
    firms <- sort(unique(df[,firm.var.name]))
    if (length(firms) <= 1) stop("Must be at least 2 firms in df.")
    
    df <- df[order(df[,time.var.name], df[,firm.var.name]),]
    
    results.df <- data.frame(time = NA, firm = NA,
                             technical.change = NA,
                             efficiency.change = NA,
                             productivity.change = NA)
    names(results.df) <- c(time.var.name, firm.var.name,
                           "technical.change", "efficiency.change", "productivity.change")
    
    #Run analysis
    #First subset first year
    df.0 <- df[df[,time.var.name] == years[1],]
    if (!all(df.0[,firm.var.name] == firms)) stop("Method currently only built for balanced panels.")
    
    X.0 <- as.matrix(df.0[,input.var.names])
    y.0 <- df.0[,output.var.name]
    
    reflected.data.0 <- reflect.data(X.0, y.0)
    X.0.eval <- reflected.data.0$X.reflected
    y.0.eval <- reflected.data.0$y.reflected
    
    #Loop through years
    for (i in 2:num.years){
      year.difference <- years[i] - years[i-1]
      year.scaling.factor <- 1 / year.difference
      
      df.1 <- df[df[,time.var.name] == years[i],]
      if (!all(df.0[,firm.var.name] == firms)) stop("Method currently only built for balanced panels.")
      
      X.1 <- as.matrix(df.1[,input.var.names])
      y.1 <- df.1[,output.var.name]
      
      reflected.data.1 <- reflect.data(X.1, y.1)
      X.1.eval <- reflected.data.1$X
      y.1.eval <- reflected.data.1$y
      
      X.constrained <- rbind(X.0, X.1)
      
      #D.s.t is efficiency of production plan in year t relative to boundary in year s
      D.0.0 <- fit.boundary(X.0.eval, y.0.eval, X.0, y.0, X.constrained, X.0, y.0, method = method)
      D.0.1 <- fit.boundary(X.0.eval, y.0.eval, X.0, y.0, X.constrained, X.1, y.1, method = method)
      D.1.0 <- fit.boundary(X.1.eval, y.1.eval, X.1, y.1, X.constrained, X.0, y.0, method = method)
      D.1.1 <- fit.boundary(X.1.eval, y.1.eval, X.1, y.1, X.constrained, X.1, y.1, method = method)
      
      #Calculate efficiency, technical changes
      efficiency.change <- D.1.1$efficiency / D.0.0$efficiency
      technical.change <- sqrt((D.0.1$efficiency / D.1.1$efficiency) * 
                                 (D.0.0$efficiency / D.1.0$efficiency))
      
      #Organize results and bind to results.df
      temp.df <- df.1[,c(time.var.name, firm.var.name)]
      temp.df$efficiency.change <- efficiency.change ^ year.scaling.factor
      temp.df$technical.change <- technical.change ^ year.scaling.factor
      temp.df$productivity.change <- (efficiency.change * technical.change) ^ year.difference
      
      results.df <- rbind(results.df, temp.df)
      
      #Set old vars to new for next round
      X.0 <- X.1
      y.0 <- y.1
      X.0.eval <- X.1.eval
      y.0.eval <- y.1.eval
    }
    results.df <- results.df[-1,]
    
    return(results.df)
  }
