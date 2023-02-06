{

//Function to get the  bootstrap uncertainty, Data stored on array data, N subsamples, mean0=the mean of the original sample 
Double_t bootstr(Double_t data[10000], Int_t N, Double_t mean0){
Double_t sum=0;
Int_t p=0;
Int_t n=10000/N; //number of entries per subsample
for (Int_t j=0;j<N;j++){
	Double_t sumj=0, meanv2j=0;
	for (Int_t s=0;s<n;s++){
		p=j*n+s;	
		sumj=sumj+data[p];
	}
	meanv2j=sumj/n;
	sum=sum+(meanv2j-mean0)*(meanv2j-mean0);
}
Double_t var=sum/(N-1);
Double_t uncertainty=sqrt(var/N);
return uncertainty;
}

//Defining function and histogram
TF1 * func = new TF1("func","(x/([1]^2))*exp(-(x^2+[0]^2)/(2*[1]^2))*ROOT::Math::cyl_bessel_i(0,x*[0]/([1]^2))",0,1);
func->SetParameter(0,0.00352672807); // I suppose a>0
func->SetParameter(1,0.0352672807);
TH1F * histo = new TH1F("histo","Sample",100,0,1);

//sampling
Int_t N=10000;
Double_t sampl[10000];
for (Int_t i=0;i<N;i++) {
Double_t v = func->GetRandom(0,1);
sampl[i]=v*v;
histo->Fill(sampl[i]);
}


cout<<histo->GetMean()<<endl; //mean of histogram: 0.00252635
cout<<histo->GetMeanError()<<endl; //error of histogram: 2.55221e-05
cout<<bootstr(sampl,10,histo->GetMean())<<endl; //bootstrap error 10 subsamples 2.57174e-05
cout<<bootstr(sampl,20,histo->GetMean())<<endl; //bootstrap error 20 subsamples 2.08996e-05 
cout<<bootstr(sampl,100,histo->GetMean())<<endl; //bootstrap error 100 subsamples 2.28325e-05
delete histo;
delete func;
}
