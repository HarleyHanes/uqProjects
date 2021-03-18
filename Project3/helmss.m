function ss=helmss(params,data)

xdata = data.xdata;
udata = data.ydata;


alpha1=params(1);
alpha11=params(2);
alpha111=params(3);

uvals_data=alpha1.*xdata.^2+alpha11.*xdata.^4+alpha111.*xdata.^6;

res = udata - uvals_data;
ss = res'*res;