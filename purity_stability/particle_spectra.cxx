#include "eicsmear/erhic/EventPythia.h"

#include "Riostream.h"
#include "TApplication.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TRint.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TLine.h"
#include "TGraph.h"
#include "TF1.h"
#include "TChain.h"
#include "TPaveText.h"

//#include "PadMxN.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

//Globals
const double Mp(0.9383);
const double Me(0.511E-3);

//---------------------
int main(int argc, char **argv){

    #ifdef WITHRINT
    TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
    #else
    TApplication *myapp = new TApplication("myapp",0,0);
    #endif

    //Eta bin for all electrons (black), positrons (green), photons (red), postive+negative pions (blue)
    int nEta = 1; //Number of eta bins
    double Eta_low[] = {-4.5};
    double Eta_hi[] =  {-3.5};

    TH1* h_elec[nEta]; //All electrons 
    TH1* h_pos[nEta];  //Positrons
    TH1* h_phot[nEta];  //Photons
    TH1* h_pi[nEta];   //Charged pions

    //Veto on Q2 tagger
    TH1* h_elec1[nEta]; //All electrons 
    TH1* h_pos1[nEta];  //Positrons
    TH1* h_phot1[nEta];  //Photons
    TH1* h_pi1[nEta];   //Charged pions

    //Veto on Q2 tagger + Loose cut of (E-pz)_tot > 18 GeV
    TH1* h_elec2[nEta]; //All electrons 
    TH1* h_pos2[nEta];  //Positrons
    TH1* h_phot2[nEta];  //Photons
    TH1* h_pi2[nEta];   //Charged pions

    //use constant log binning
    double p_min = 1E-1;
    double p_max = 40;
	const int nbins = 100;
	double log_bw = (log10(p_max) - log10(p_min))/((double)nbins);
    double log_div;
    double p_bins[nbins+1];
    for(int i=0;i<nbins+1;i++){
		log_div = log10(p_min) + (i*log_bw);
		p_bins[i] = pow(10,log_div);
	}

    for(int ihist=0;ihist<nEta;ihist++){

        h_elec[ihist] = new TH1D(Form("h_elec[%d]",ihist),"",nbins,p_bins);
        h_elec[ihist]->SetLineColor(kBlack);h_elec[ihist]->SetLineWidth(2);

        h_pos[ihist] = new TH1D(Form("h_pos[%d]",ihist),"",nbins,p_bins);
        h_pos[ihist]->SetLineColor(kGreen);h_pos[ihist]->SetLineWidth(2);

        h_phot[ihist] = new TH1D(Form("h_phot[%d]",ihist),"",nbins,p_bins);
        h_phot[ihist]->SetLineColor(kRed);h_phot[ihist]->SetLineWidth(2);

        h_pi[ihist] = new TH1D(Form("h_pi[%d]",ihist),"",nbins,p_bins);
        h_pi[ihist]->SetLineColor(kBlue);h_pi[ihist]->SetLineWidth(2);


        h_elec1[ihist] = new TH1D(Form("h_elec1[%d]",ihist),"",nbins,p_bins);
        h_elec1[ihist]->SetLineColor(kBlack);h_elec1[ihist]->SetLineWidth(2);

        h_pos1[ihist] = new TH1D(Form("h_pos1[%d]",ihist),"",nbins,p_bins);
        h_pos1[ihist]->SetLineColor(kGreen);h_pos1[ihist]->SetLineWidth(2);

        h_phot1[ihist] = new TH1D(Form("h_phot1[%d]",ihist),"",nbins,p_bins);
        h_phot1[ihist]->SetLineColor(kRed);h_phot1[ihist]->SetLineWidth(2);

        h_pi1[ihist] = new TH1D(Form("h_pi1[%d]",ihist),"",nbins,p_bins);
        h_pi1[ihist]->SetLineColor(kBlue);h_pi1[ihist]->SetLineWidth(2);


        h_elec2[ihist] = new TH1D(Form("h_elec2[%d]",ihist),"",nbins,p_bins);
        h_elec2[ihist]->SetLineColor(kBlack);h_elec2[ihist]->SetLineWidth(2);

        h_pos2[ihist] = new TH1D(Form("h_pos2[%d]",ihist),"",nbins,p_bins);
        h_pos2[ihist]->SetLineColor(kGreen);h_pos2[ihist]->SetLineWidth(2);

        h_phot2[ihist] = new TH1D(Form("h_phot2[%d]",ihist),"",nbins,p_bins);
        h_phot2[ihist]->SetLineColor(kRed);h_phot2[ihist]->SetLineWidth(2);

        h_pi2[ihist] = new TH1D(Form("h_pi2[%d]",ihist),"",nbins,p_bins);
        h_pi2[ihist]->SetLineColor(kBlue);h_pi2[ihist]->SetLineWidth(2);
    }

    //Load ROOT Files
    erhic::EventPythia *event(NULL); //Event Class
    erhic::ParticleMC *particle(NULL); //Particle Class

    TChain *tree = new TChain("EICTree");

    //Using files created with my installation
    for(int i=0;i<300;i++){
        tree->Add(Form("/gpfs02/eic/baraks/pythia/outfiles/other_studies/18_275/fullQ2/ep_minbias_%d.root",i));
    }

    tree->SetBranchAddress("event",&event);

    //Particle variables
    int nParticles(0);
    int id(0);
    int status(0);
    int orig(0);

    //Calculate Generated Luminosity
    int nevents = tree->GetEntries();
    double cross_tot = 181.8E9; //Total Cross Section in fb
    double lum = ( (double) nevents)/cross_tot; //Luminosity in fb^-1

    cout<<"-------------------------------"<<endl;
    cout<<"PYTHIA6 Simulation:"<<endl;
    cout<<"Total Number of Events = "<<nevents<<endl;
    cout<<"Integrated Luminosity = "<<lum<<" fb^-1"<<endl<<endl;

    //Loop over events
    for(int iEvent=0;iEvent<nevents;iEvent++){
        
        if(iEvent%10000==0) cout<<"Events Analysed = "<<iEvent<<"!"<<endl;
        tree->GetEntry(iEvent);

        //Get Number of tracks
        nParticles = event->GetNTracks();

        //Cut variables
        double Q2_veto = false; //Set to true if an electron is found in low Q2 tagger
        double empztot = 0;

        // First loop over all particles in event
        for(int iParticle=0;iParticle<nParticles;iParticle++){

            particle = event->GetTrack(iParticle);

            id = (int) particle->Id();
            status = (int) particle->GetStatus();
            orig = (int) particle->GetParentIndex();

            //Veto on low Q2 tagger (based on figure 11.120b)     
            auto Q2_true = event->GetTrueQ2();
            //auto x_true = event->GetTrueX();
            auto y_true = event->GetTrueY();

            if( log10(y_true)<-1.1 && log10(Q2_true)>-2. && log10(Q2_true)<-1.25 ) //1st region
                Q2_veto = true;
            else if( log10(y_true)>-1.1 && log10(y_true)<-0.3 && log10(Q2_true)>-7.5 && log10(Q2_true)<-1.25 ) //2nd region
                Q2_veto = true;

            //Total E-pz for event
            if( status==1 && particle->GetEta()>-4.5 && particle->GetEta()<4.0 ){
                empztot+=( particle->GetE()-particle->GetPz() );
            }
        }

        // Loop again over all particles in event
        for(int iParticle=0;iParticle<nParticles;iParticle++){
      
            particle = event->GetTrack(iParticle);

            id = (int) particle->Id();
            status = (int) particle->GetStatus();
            orig = (int) particle->GetParentIndex();

            if(id==11 && status==1 && orig==3){ //Scattered Electrons
                for(int iEta=0;iEta<nEta;iEta++){
                    if(particle->GetEta()>Eta_low[iEta] && particle->GetEta()<Eta_hi[iEta]){
                            h_elec[iEta]->Fill(particle->GetP());
                            if(!Q2_veto) h_elec1[iEta]->Fill(particle->GetP());
                            if(!Q2_veto && empztot>18) h_elec2[iEta]->Fill(particle->GetP());
                    }
                }
            }

            if(id==11 && status==1 && orig!=3){ //Non-scattered Electrons
                for(int iEta=0;iEta<nEta;iEta++){
                    if(particle->GetEta()>Eta_low[iEta] && particle->GetEta()<Eta_hi[iEta]){
                        h_elec[iEta]->Fill(particle->GetP()); //Fill same histogram as scattered electron
                        if(!Q2_veto) h_elec1[iEta]->Fill(particle->GetP());
                        if(!Q2_veto && empztot>18) h_elec2[iEta]->Fill(particle->GetP());
                    }
                }
            }

            if(id==-11 && status==1){ //Positrons
                for(int iEta=0;iEta<nEta;iEta++){
                    if(particle->GetEta()>Eta_low[iEta] && particle->GetEta()<Eta_hi[iEta]){
                        h_pos[iEta]->Fill(particle->GetP());
                        if(!Q2_veto) h_pos1[iEta]->Fill(particle->GetP());
                        if(!Q2_veto && empztot>18) h_pos2[iEta]->Fill(particle->GetP());
                    }
                }
            }

            if(id==-211 && status==1){ //Negative Pions
                for(int iEta=0;iEta<nEta;iEta++){
                    if(particle->GetEta()>Eta_low[iEta] && particle->GetEta()<Eta_hi[iEta]){
                        h_pi[iEta]->Fill(particle->GetP());
                        if(!Q2_veto) h_pi1[iEta]->Fill(particle->GetP());
                        if(!Q2_veto && empztot>18) h_pi2[iEta]->Fill(particle->GetP());
                    }
                }
            }

            if(id==211 && status==1){ //Positive Pions
                for(int iEta=0;iEta<nEta;iEta++){
                    if(particle->GetEta()>Eta_low[iEta] && particle->GetEta()<Eta_hi[iEta]){
                        h_pi[iEta]->Fill(particle->GetP()); //Fill same histogram as negative pions
                        if(!Q2_veto) h_pi1[iEta]->Fill(particle->GetP());
                        if(!Q2_veto && empztot>18) h_pi2[iEta]->Fill(particle->GetP());
                    }
                }
            }

            if(id==22 && status==1){ //Photons
                for(int iEta=0;iEta<nEta;iEta++){
                    if(particle->GetEta()>Eta_low[iEta] && particle->GetEta()<Eta_hi[iEta]){
                        h_phot[iEta]->Fill(particle->GetP());
                        if(!Q2_veto) h_phot1[iEta]->Fill(particle->GetP());
                        if(!Q2_veto && empztot>18) h_phot2[iEta]->Fill(particle->GetP());
                    }
                }
            }

        } //closes particle loop
    } //closes event loop

    //Scale histograms
    for(int ihist=0;ihist<nEta;ihist++){
        h_elec[ihist]->Scale(1./lum);
        h_pos[ihist]->Scale(1./lum);
        h_pi[ihist]->Scale(1./lum);
        h_phot[ihist]->Scale(1./lum);

        h_elec1[ihist]->Scale(1./lum);
        h_pos1[ihist]->Scale(1./lum);
        h_pi1[ihist]->Scale(1./lum);
        h_phot1[ihist]->Scale(1./lum);

        h_elec2[ihist]->Scale(1./lum);
        h_pos2[ihist]->Scale(1./lum);
        h_pi2[ihist]->Scale(1./lum);
        h_phot2[ihist]->Scale(1./lum);
        
    }

    //Draw plots
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    gStyle->SetLabelSize(0.04,"X");
    gStyle->SetLabelSize(0.04,"Y");
    //gStyle->SetLabelOffset(0.01,"X");
    //gStyle->SetLabelOffset(0.01,"Y");
    gStyle->SetTitleXSize(0.055);
    gStyle->SetTitleXOffset(0.85);
    gStyle->SetTitleYSize(0.055);
    gStyle->SetTitleYOffset(0.85);

    //Plot 1
    TCanvas *c1 = new TCanvas("c1");
	TH2 *hframe_1 = new TH2F("hframe_1","",100,p_min,p_max,100,0.2,1e12);

    gPad->SetLogx();gPad->SetLogy();
    gPad->SetTickx();gPad->SetTicky();

    hframe_1->Draw();
    hframe_1->GetXaxis()->SetTitle("Momentum [GeV/c]");hframe_1->GetYaxis()->SetTitle("Particles / bin / fb^{-1}");
    hframe_1->GetXaxis()->CenterTitle(1);hframe_1->GetYaxis()->CenterTitle(1);
	//hframe_1->GetXaxis()->SetLabelFont(63);hframe_1->GetYaxis()->SetLabelFont(63);
    //hframe_1->GetXaxis()->SetLabelSize(25);hframe_1->GetYaxis()->SetLabelSize(25);
    //hframe_1->GetXaxis()->SetLabelOffset(0.01);hframe_1->GetYaxis()->SetLabelOffset(0.01);
    //hframe_1->GetXaxis()->SetTitleSize(40);hframe_1->GetXaxis()->SetTitleOffset(2.5); 
    //hframe_1->GetYaxis()->SetTitleSize(40);hframe_1->GetYaxis()->SetTitleOffset(3.0);
        
    h_elec[0]->Draw("hist same");
    h_pos[0]->Draw("hist same");
    h_phot[0]->Draw("hist same");
    h_pi[0]->Draw("hist same");

    //Draw Text
    TPaveText* pave1_1 = new TPaveText(0.25,0.8,0.35,0.9,"NDCNB");
    pave1_1->AddText("18 GeV e on 275 GeV p"); //Hard code the beam energy for now
	pave1_1->SetFillStyle(4000);
    pave1_1->SetTextFont(63);pave1_1->SetTextSize(25);
    pave1_1->Draw();

    TPaveText* pave1_2 = new TPaveText(0.65,0.85,0.8,0.9,"NDCNB");
    pave1_2->AddText("All Electrons");
	pave1_2->SetFillStyle(4000);
    pave1_2->SetTextFont(63);pave1_2->SetTextSize(18);
    pave1_2->SetTextColor(kBlack);
    pave1_2->Draw();

    TPaveText* pave1_3 = new TPaveText(0.65,0.8,0.8,0.85,"NDCNB");
    pave1_3->AddText("Charged Pions");
	pave1_3->SetFillStyle(4000);
    pave1_3->SetTextFont(63);pave1_3->SetTextSize(18);
    pave1_3->SetTextColor(kBlue);
    pave1_3->Draw();

    TPaveText* pave1_4 = new TPaveText(0.65,0.77,0.8,0.8,"NDCNB");
    pave1_4->AddText("Positrons");
	pave1_4->SetFillStyle(4000);
    pave1_4->SetTextFont(63);pave1_4->SetTextSize(18);
    pave1_4->SetTextColor(kGreen);
    pave1_4->Draw();

    TPaveText* pave1_5 = new TPaveText(0.65,0.72,0.8,0.75,"NDCNB");
    pave1_5->AddText("Photons");
	pave1_5->SetFillStyle(4000);
    pave1_5->SetTextFont(63);pave1_5->SetTextSize(18);
    pave1_5->SetTextColor(kRed);
    pave1_5->Draw();

    TPaveText *pave2[nEta];
    double p2_xl[] = {0.2};
    double p2_xh[] = {0.4};
    double p2_yl[] = {0.75};
    double p2_yh[] = {0.8};

    for(int ipave=0;ipave<nEta;ipave++){
        pave2[ipave] = new TPaveText(p2_xl[ipave],p2_yl[ipave],p2_xh[ipave],p2_yh[ipave],"NDCNB");
        pave2[ipave]->AddText(Form("%.1f < #eta < %.1f",Eta_low[ipave],Eta_hi[ipave]));
	    pave2[ipave]->SetFillStyle(4000);
        pave2[ipave]->SetTextFont(63);pave2[ipave]->SetTextSize(25);
        pave2[ipave]->Draw();
    }

    //Create vertical line at minimum momentum values --
    //defined as electron minimum electron momentum in that eta range satifying y<0.95
    double mom_low[] = {0.9};
    TLine *linea[nEta];

    TPaveText *pave3[nEta];
    double p3_xl[] = {0.4};
    double p3_xh[] = {0.45};
    double p3_yl[] = {0.7};
    double p3_yh[] = {0.75};

    for(int iline=0;iline<nEta;iline++){
        linea[iline] = new TLine(mom_low[iline],0.2,mom_low[iline],h_elec[0]->GetMaximum());
        linea[iline]->SetLineColor(kOrange);
        linea[iline]->SetLineWidth(2);
        linea[iline]->SetLineStyle(2);
        linea[iline]->Draw();

        pave3[iline] = new TPaveText(p3_xl[iline],p3_yl[iline],p3_xh[iline],p3_yh[iline],"NDCNB");
        pave3[iline]->AddText("P_{min.}");
	    pave3[iline]->SetFillStyle(4000);
        pave3[iline]->SetTextFont(63);pave3[iline]->SetTextSize(15);
        pave3[iline]->SetTextColor(kOrange);
        pave3[iline]->Draw();
    }

    c1->Modified();c1->Update();

    //Plot 2
    TCanvas *c2 = new TCanvas("c2");
    gPad->SetLogx();gPad->SetLogy();
    gPad->SetTickx();gPad->SetTicky();
    hframe_1->Draw();

    h_elec1[0]->Draw("hist same");
    h_pos1[0]->Draw("hist same");
    h_phot1[0]->Draw("hist same");
    h_pi1[0]->Draw("hist same");

    pave1_1->Draw();
    pave1_2->Draw();
    pave1_3->Draw();
    pave1_4->Draw(); 
    pave1_5->Draw();
    pave2[0]->Draw();
    linea[0]->Draw();pave3[0]->Draw();

    c2->Modified();c2->Update();

    //Plot 3
    TCanvas *c3 = new TCanvas("c3");
    gPad->SetLogx();gPad->SetLogy();
    gPad->SetTickx();gPad->SetTicky();
    hframe_1->Draw();

    h_elec2[0]->Draw("hist same");
    h_pos2[0]->Draw("hist same");
    h_phot2[0]->Draw("hist same");
    h_pi2[0]->Draw("hist same");

    pave1_1->Draw();
    pave1_2->Draw();
    pave1_3->Draw();
    pave1_4->Draw(); 
    pave1_5->Draw();
    pave2[0]->Draw();
    linea[0]->Draw();pave3[0]->Draw();

    c3->Modified();c3->Update();

    //Print to file
    gROOT->ProcessLine("c1->Print(\"plots/particle_spectra.pdf[\");");
    gROOT->ProcessLine("c1->Print(\"plots/particle_spectra.pdf\");");
    gROOT->ProcessLine("c2->Print(\"plots/particle_spectra.pdf\");");
    gROOT->ProcessLine("c3->Print(\"plots/particle_spectra.pdf\");");
    gROOT->ProcessLine("c3->Print(\"plots/particle_spectra.pdf]\");");

    //Write histograms to file
    /* TList *l = new TList();
    for(int ihist=0;ihist<nEta;ihist++){
        l->Add(h_elec[ihist]);
        l->Add(h_pos[ihist]);
        l->Add(h_pi[ihist]);
        l->Add(h_phot[ihist]);

        l->Add(h_elec1[ihist]);
        l->Add(h_pos1[ihist]);
        l->Add(h_pi1[ihist]);
        l->Add(h_phot1[ihist]);

        l->Add(h_elec2[ihist]);
        l->Add(h_pos2[ihist]);
        l->Add(h_pi2[ihist]);
        l->Add(h_phot2[ihist]);
    }

    l->Write("histlist", TObject::kSingleKey); */

    //ROOT file to save histograms
    TFile *fout = new TFile("out_hists.root","RECREATE");

    TH1* elec = (TH1*) h_elec[0]->Clone("elec");
    TH1* pos = (TH1*) h_pos[0]->Clone("pos");
    TH1* pi = (TH1*) h_pi[0]->Clone("pi");
    TH1* phot = (TH1*) h_phot[0]->Clone("phot");

    TH1* elec1 = (TH1*) h_elec1[0]->Clone("elec1");
    TH1* pos1 = (TH1*) h_pos1[0]->Clone("pos1");
    TH1* pi1 = (TH1*) h_pi1[0]->Clone("pi1");
    TH1* phot1 = (TH1*) h_phot1[0]->Clone("phot1");

    TH1* elec2 = (TH1*) h_elec2[0]->Clone("elec2");
    TH1* pos2 = (TH1*) h_pos2[0]->Clone("pos2");
    TH1* pi2 = (TH1*) h_pi2[0]->Clone("pi2");
    TH1* phot2 = (TH1*) h_phot2[0]->Clone("phot2");

    elec->Write();pos->Write();pi->Write();phot->Write();
    elec1->Write();pos1->Write();pi1->Write();phot1->Write();
    elec2->Write();pos2->Write();pi2->Write();phot2->Write();

    fout->Close();

    myapp->Run();
    return 0;
}
