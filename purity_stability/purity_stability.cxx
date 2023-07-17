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
#include "TRandom3.h"

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
 
    //Q2 Binning                                                                        
    double Q2_min = 1E-2;
    double Q2_max = 1E1;
    const int nbins_Q2 = 20;
    double log_bw_Q2 = (log10(Q2_max) - log10(Q2_min))/(nbins_Q2); //Determine bin width                                                         
    double log_Q2_div;
    double Q2_bins[nbins_Q2+1];
    for(int i=0;i<nbins_Q2+1;i++){
        log_Q2_div = log10(Q2_min) + (i*log_bw_Q2);
        Q2_bins[i] = pow(10,log_Q2_div);
    }

    //x Binning                                                                                                                       
    double x_min = 1E-6;
    double x_max = 1;
    const int nbins_x = 30;
    double log_bw_x = (log10(x_max) - log10(x_min))/(nbins_x); //Determine bin width                                                             
    double log_x_div;
    double x_bins[nbins_x+1];
    for(int i=0;i<nbins_x+1;i++){
        log_x_div = log10(x_min) + (i*log_bw_x);
        x_bins[i] = pow(10,log_x_div);
    }

    //Index holding variables
    int genbin(0), recbin(0);

    //Random number generator
    TRandom3 *myrand = new TRandom3(0);

    //Output ROOT File
    TFile *fout = new TFile("purity_stability_hists.root","RECREATE");

    //Histograms
    TH2D *h0 = new TH2D("h0","Generated",nbins_x,x_bins,nbins_Q2,Q2_bins);
    h0->GetXaxis()->SetTitle("x");h0->GetXaxis()->CenterTitle();
    h0->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h0->GetYaxis()->CenterTitle();

    TH2D *h1 = new TH2D("h1","Reconstructed",nbins_x,x_bins,nbins_Q2,Q2_bins);
    h1->GetXaxis()->SetTitle("x");h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h1->GetYaxis()->CenterTitle();

    TH2D *h1_1 = new TH2D("h1_1","Electron Reconstruction Purity",nbins_x,x_bins,nbins_Q2,Q2_bins);
    h1_1->GetXaxis()->SetTitle("x");h1_1->GetXaxis()->CenterTitle();
    h1_1->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h1_1->GetYaxis()->CenterTitle();
    
    TH2D *h1_2 = new TH2D("h1_2","Electron Reconstruction Stability",nbins_x,x_bins,nbins_Q2,Q2_bins);
    h1_2->GetXaxis()->SetTitle("x");h1_2->GetXaxis()->CenterTitle();
    h1_2->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");h1_2->GetYaxis()->CenterTitle();
    
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
    int nevents = 1E7;//tree->GetEntries();
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

        // Loop over all particles in event
        for(int iParticle=0;iParticle<nParticles;iParticle++){

            particle = event->GetTrack(iParticle);

            id = (int) particle->Id();
            status = (int) particle->GetStatus();
            orig = (int) particle->GetParentIndex();
              
            auto Q2_true = event->GetTrueQ2();
            auto x_true = event->GetTrueX();
            //auto y_true = event->GetTrueY();

            //Scattered electron in fdc acceptance
            if( id==11 && status==1 && orig==3 && particle->GetEta()>-4.6 && particle->GetEta()<-3.6 ){

                //Fill generated histograms
                genbin = h0->Fill(x_true,Q2_true);
                
                //Smear electron particle energy and angles
                auto E_true = particle->GetE();
                auto E_sigma = E_true*sqrt( pow(0.17/sqrt(E_true),2) + pow(0.02,2) );
                auto E_rec = myrand->Gaus(E_true,E_sigma);

                auto rad_true = tan(particle->GetTheta()) * -3070.; //in mm
                auto xpos_true = rad_true * cos(particle->GetPhi());
                auto ypos_true = rad_true * sin(particle->GetPhi());

                auto pos_sigma = 3./sqrt(E_true); //in mm

                auto xpos_rec = myrand->Gaus(xpos_true,pos_sigma);
                auto ypos_rec = myrand->Gaus(ypos_true,pos_sigma);
                auto rad_rec = sqrt( pow(xpos_rec,2) + pow(ypos_rec,2) );
                auto theta_rec = atan2(rad_rec,-3070.);

                if(E_rec>0){
                    //Calculate reconstructed x and Q2
                    auto Q2_rec = 4. * 18. * E_rec * cos(theta_rec/2.) * cos(theta_rec/2.);
                    auto y_rec = 1. - ( (E_rec/(2.*18.)) * (1. - cos(theta_rec)) );
                    auto x_rec = Q2_rec / (4.*18.*275.*y_rec);

                    //Fill reconstructed histogram
                    recbin = h1->Fill(x_rec,Q2_rec);

                    //Purity/Stability histograms
                    if(genbin==recbin){
                        h1_1->Fill(x_rec,Q2_rec);
                        h1_2->Fill(x_rec,Q2_rec);
                    }
                }

                break; //Only need scattered electron
            }
        } //closes particle loop
    } //closes event loop

    //Divide histograms
    h1_1->Divide(h1);
    h1_2->Divide(h0);

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

    //Define y=constant functions
    double y_max = 0.95;double y_min = 1e-2;double s_cm = 4.*18.*275.;

    TF1 *f_ymin = new TF1("f_ymin","x*[0]*[1]",1e-6,1);
    f_ymin->SetLineColor(kRed);f_ymin->SetLineStyle(2);f_ymin->SetLineWidth(2);
    f_ymin->SetParameters(s_cm,y_min);

    TF1 *f_ymax = new TF1("f_ymax","x*[0]*[1]",1e-6,1);
    f_ymax->SetLineColor(kRed);f_ymax->SetLineStyle(2);f_ymax->SetLineWidth(2);
    f_ymax->SetParameters(s_cm,y_max);

    //Make Latex
    TLatex *tex_energy = new TLatex();
    TLatex *tex_ymin = new TLatex();
    TLatex *tex_ymax = new TLatex();

    tex_energy->SetText(1.5e-6,8,"18 GeV e on 275 GeV p");
    tex_energy->SetTextColor(kBlack);
    tex_energy->SetTextSize(0.03);

    tex_ymin->SetText(2e-2,2,Form("y = %.2f",y_min));
    tex_ymin->SetTextColor(kRed);
    tex_ymin->SetTextSize(0.035);tex_ymin->SetTextAngle(70);

    tex_ymax->SetText(2e-4,2,Form("y = %.2f",y_max));
    tex_ymax->SetTextColor(kRed);
    tex_ymax->SetTextSize(0.035);tex_ymax->SetTextAngle(70);
    
    //Make Plot
    TCanvas *c1 = new TCanvas("c1");
    
    c1->Divide(2,1);
    c1->cd(1);gPad->SetTopMargin(0.12);gPad->SetBottomMargin(0.12);gPad->SetRightMargin(0.12);gPad->SetLeftMargin(0.12);
    gPad->SetLogx();gPad->SetLogy();h1_1->Draw("colz");
    f_ymin->Draw("same");f_ymax->Draw("same");tex_ymin->Draw();tex_ymax->Draw();
    
    c1->cd(2);gPad->SetTopMargin(0.12);gPad->SetBottomMargin(0.12);gPad->SetRightMargin(0.12);gPad->SetLeftMargin(0.12);
    gPad->SetLogx();gPad->SetLogy();h1_2->Draw("colz");
    f_ymin->Draw("same");f_ymax->Draw("same");tex_ymin->Draw();tex_ymax->Draw();
    tex_energy->Draw();

    //Print to file
    c1->Print("plots/purity_stability.pdf");

    //Write output ROOT file
    h1_1->Write();h1_2->Write();
    fout->Close();
    
    myapp->Run();
    return 0;
}
