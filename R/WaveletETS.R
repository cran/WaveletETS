
#' @title Wavelet Based Error Trend Seasonality Model
#' @param ts Time Series Data
#' @param split_ratio Training and Testing Split
#' @param wlevels Number of Wavelet Levels
#' @import dplyr Metrics tseries stats wavelets forecast caretForecast
#' @return
#' \itemize{
#'   \item Train_actual: Actual train series
#'   \item Test_actual: Actual test series
#'   \item Train_fitted: Fitted train series
#'   \item Test_predicted: Predicted test series
#'   \item Accuracy: RMSE and MAPE of the model
#' }
#'
#' @export
#'
#' @examples
#' library("WaveletETS")
#' data<- rnorm(100,100, 10)
#' WG<-WaveletETS(ts=data)
#' @references
#' \itemize{
#'\item Aminghafari, M. and Poggi, J.M. 2012. Nonstationary time series forecasting using wavelets and kernel smoothing. Communications in Statistics-Theory and Methods, 41(3),485-499.

#' \item Paul, R.K. A and Anjoy, P. 2018. Modeling fractionally integrated maximum temperature series in India in presence of structural break. Theory and Applied Climatology 134, 241–249.

#' }
WaveletETS<-function(ts,split_ratio=0.8,wlevels=3){
  ntest<-round(length(ts)*(1-split_ratio), digits = 0)
  Split1 <- caretForecast::split_ts(as.ts(ts), test_size = ntest)
  train_data1 <- Split1$train
  test_data1 <- Split1$test
  Wvlevels<-wlevels
  mraout <- wavelets::modwt(as.vector(ts), filter="haar", n.levels=Wvlevels)
  WaveletSeries <- cbind(do.call(cbind,mraout@W),mraout@V[[Wvlevels]])
  ts_fitted<-NULL
  ts_foreast<-NULL

  TSModel<-function(ts,split_ratio,method){
    ntest<-round(length(ts)*(1-split_ratio), digits = 0)
    Split <- caretForecast::split_ts(as.ts(ts), test_size = ntest)
    train_data <- Split$train
    test_data <- Split$test
    model<- ets(train_data)
    ts_fitted<-as.vector(model$fitted)
    ts_foreast<-as.vector(forecast(model, h=ntest)$mean)
    return(list(Model=model,Train_actual=train_data,Test_actual=test_data, Train_fitted=ts_fitted,Test_predicted=ts_foreast))
  }

  for (j in 1:ncol(WaveletSeries)) {
    w<-as.ts(WaveletSeries[,j])
    model<-TSModel(ts=w,split_ratio=split_ratio,method="gbm")
    ts_fitted<-cbind(ts_fitted,model$Train_fitted)
    ts_foreast<-cbind(ts_foreast,model$Test_predicted)
  }

  trainf <- apply(ts_fitted,1,sum)
  testf <- apply(ts_foreast,1,sum)

  RMSE<-c(Train=Metrics::rmse(train_data1,trainf),Test=Metrics::rmse(test_data1,testf))
  MAPE<-c(Train=Metrics::mape(train_data1,trainf),Test=Metrics::mape(test_data1,testf))
  accuracy<-rbind(RMSE,MAPE)
  return(list(Train_actual=train_data1,Test_actual=test_data1,Train_fitted=trainf,Test_predicted=testf, Accuracy=accuracy))
}

