#include <iostream>
using namespace std;




//Simple macro for 
//1. estimate of pile-up within one bunch crossing with Poissonian distribution in ALICE
//2. estimate of several bunch crossings leading to collision within integration time window of a detector (by default drift time of ALICE TPC)

//Author: Michael Winn, 15.3.2013

///Assumptions:
//For 2. "equal-distance" distribution of bunches over ring using the number of colliding bunches and the revolution frequency is assumed 
//(for integration times in the order of one orbit (order of 100 mikroseconds), right assumption, but one has to take into account that the bunch crossings are not homogeneously distributed 
//over the integration time, for smaller integration times, this is only a time average over one orbit), 
//or the bunch spacing directly (this assumption only reasonable approximation, when machine pretty full and nearly all  bunches colliding in ALICE for large integration time
// or when one is interested in small integration times (small w.r.t. 1/revolution_frequency = 89 mikroseconds or better ) respectively the number of pile-up within the middle of one bunch train ), so the the latter version can serve as an upper limit of the effect of pile-up

//same intensity for all colliding bunches assumed

//ATTENTION: interaction frequency might change considerably within one fill or even within one ALICE-run!

//ATTENTION: trigger can bias strongly the occurence of pile-up in recorded events (e.g. high-mult. trigger from SPD might be strongly biased towards pile-up)

//web-pages to be considered for validity of assumptions:
//http://lpc.web.cern.ch/lpc/fillingschemes.htm
//https://alice-logbook.cern.ch/logbook/date_online.php?p_cont=fd&p_fill=3474&p_tab=fls
//https://alice-logbook.cern.ch/logbook/date_online.php?p_cont=rund&p_run=195351&p_tab=gi
//http://alimonitor.cern.ch/configuration/index.jsp?partition=LHC13b&pass=1&raw_run=&filling_scheme=&filling_config=&fillno=&energy=&intensity_per_bunch=&mu=&interacting_bunches=&noninteracting_bunches_beam_1=&noninteracting_bunches_beam_2=&interaction_trigger=&rate=&beam_empty_trigger=&empty_empty_trigger=&muon_trigger=&high_multiplicity_trigger=&emcal_trigger=&calibration_trigger=&quality=&muon_quality=&comment=&field=&det_aco=&det_emc=&det_fmd=&det_hlt=&det_hmp=&det_mch=&det_mtr=&det_phs=&det_pmd=&det_spd=&det_sdd=&det_ssd=&det_tof=&det_tpc=&det_trd=&det_t00=&det_v00=&det_zdc=&hlt_mode=&changedon=





void Pile_up_Estimator(
		       Double_t interaction_rate = 10 /*in kHz*/, 
		       Int_t n_colliding_bunch_pairs = 202 /*in ALICE, 3rd number in LHC filling scheme*/,
		       Double_t bunch_spacing = 200 /*in ns, first number in filling scheme of LHC*/, 
		       Double_t integration_time = 94 /* in mikroseconds, by default drift time of TPC according to NIM A 622 (2010) 316-367, integration time of SPD 0.2 mikroseconds(A. Mastroserio and M. Ivanov), 2.6 mikroseconds for TRD (rough estimate of drift time))*/, 
		       Int_t mode = 1 /*1: using number of colliding bunches, 0: using bunch spacing*/
		       ){


  //revolution frequency of LHC in kHz, derived from LHC web page
  const Double_t revolution_frequency = 11.245;

  //estimate of mue
  Double_t mue = interaction_rate/(revolution_frequency* n_colliding_bunch_pairs);

  cout <<"estimate for mue value: " << mue << endl;

  //probability for more than 0 interactions in 1 bunch crossing
  Double_t interaction_single_bunch_crossing = 1 - exp(-mue);
  cout <<"probability for more than 0 interactions in 1 bunch crossing: " << interaction_single_bunch_crossing << endl;

  //probability for more than 1 interaction in 1 bunch crossing
  Double_t pile_up_single_bunch_crossing =  interaction_single_bunch_crossing - mue *  exp(-mue);
  cout <<"probability for more than 1 interactions in 1 bunch crossing: " <<pile_up_single_bunch_crossing << endl;
  //

  //ratio of probability with more than one collision and with thos woth one collision 
  Double_t ratio = pile_up_single_bunch_crossing /(  mue *  exp(-mue) ); 

  cout <<"ratio of probability with more than one collision and with those with one collision: " << ratio << endl;

  cout <<"fraction of events with pile-up within same bunch-crossing (assuming no offline rejection and no trigger bias): " << ratio/(1.0+ratio) << endl;

  //number of bunch crossings within considered integration time
  Int_t n_bunch_crossings = 0;

  if(mode==1){
    //appropriate formulae for long integration times and for averages
    n_bunch_crossings = (Int_t) integration_time*exp(log(10)*(-6))/(1/(revolution_frequency*1000* n_colliding_bunch_pairs) ) ;
  }else if(mode == 0){
    //appropriate  for maximal value in case of  small integration times << 1/revolution_frequency = 89 mikroseconds
    n_bunch_crossings= (Int_t)integration_time*exp(log(10)*(-6))/(bunch_spacing*exp(log(10)*(-9)) );
  
  }else{cout <<"mode not forseen!" << endl; return;}

  cout <<"number of bunch crossings within considered integration time: " << n_bunch_crossings << endl;

  //calculating probability for more than one bunch crossing each leading to at least 1 collision within the integration time
  //(using binomial distribution for more than 1 sucesses within n_bunch_crossings trials)

  Double_t probability_pile_up_different_bunches = 1 - (n_bunch_crossings)* (1- exp(-mue))* pow(exp(-mue),n_bunch_crossings-1 )
    -  pow(exp(-mue),n_bunch_crossings);


  cout << "probability for more than one bunch crossing each leading to at least 1 collision within the integration time: " << probability_pile_up_different_bunches << endl;

  //calculating probability for exactly one bunch crossing leading to at least 1 collision within the integration time
  //(using binomial distribution for more than 1 sucesses within n_bunch_crossings trials)
 
  Double_t probability_no_pile_up_different_bunches = (n_bunch_crossings)* (1- exp(-mue))* pow(exp(-mue),n_bunch_crossings-1 );

  cout << "probability for exactly one bunch crossing leading to at least 1 collision within the integration time: " << probability_no_pile_up_different_bunches<< endl;

  //ratio between probability for more than one bunch crossing each leading to at least 1 collision within the integration time and probability for exactly one bunch crossing leading to at least 1 collision within the integration time

  Double_t ratio2 = probability_pile_up_different_bunches / probability_no_pile_up_different_bunches ;
  cout << "ratio between probability for more than one bunch crossing each leading to at least 1 collision within the integration time and probability for exactly one bunch crossing leading to at least 1 collision within the integration time: " << ratio2 << endl;

  cout <<"fraction of events with pile-up from different bunch crossings (assuming no offline rejection and no trigger bias): " << ratio2/(1.0+ratio2) << endl;

  return;

}
