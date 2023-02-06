{
Double_t pi=TMath::Pi();
THStack *hs = new THStack("hs","Stacked");

//uniform acceptance
TRandom rnds1;
TH1F *h1= new TH1F("h1","",100,0.,2*pi);
for (int i=0;i<1000;i++){ // #events
	Double_t reactiveplane1=rnds1.Uniform(2*pi);
	for (int j=0;j<1000;j++){ // #particles
		h1->Fill(rnds1.Uniform(0,2*pi));
	}
}
h1->SetMinimum(0);
h1->SetMarkerColor(kBlue);
hs->Add(h1);

// non-uniform acceptance
TF1 * fvarphi = new TF1("fvarphi","1/2*pi*(1+2*0.05*cos(x-[0])+2*0.06*cos(2*(x-[0]))+2*0.07*cos(3*(x-[0]))+2*0.08*cos(4*(x-[0]))+2*0.09*cos(5*(x-[0]))+2*0.10*cos(6*(x-[0])))",0,2*pi);
TRandom rnds2;
TH1F *h2= new TH1F("h2","",100,0.,2*pi);
for (int l=0;l<1000;l++){ // #events
	fvarphi->SetParameter(0,rnds2.Uniform(2*pi));
	for (int m=0;m<1000;m++){ // #particles
		h2->Fill(fvarphi->GetRandom(0,2*pi));
	}
}
h2->SetMinimum(0);
h2->SetMarkerColor(kRed);
hs->Add(h2);
hs->Draw();
}
