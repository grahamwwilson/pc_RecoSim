//       #include "histset2enums.h" 
// bookeeping enumeration: if we do this we don't need to worry about hist pointer copies and merging
       enum th1d_ids{
		id_ptHist,
	 	id_pzHist,
		id_numpcHist,
		id_numpccutHist,
                id_ptCutHist,
		id_pzCutHist,
		id_numHGNPCHist,
                id_ptHCutHist,
                id_pzHCutHist,
	        id_rerrHGNHist,	
		id_PVndof,
		id_r25RHist,
		id_r25CHist,
		id_r25HHist,
		id_r25Hist_b2p5,
		id_r25Hist_b2p5_nowt,
		id_r25coarse,
		id_npv,
		id_pvz,
		id_pvdz,
		id_pvtrksum,
		id_sumpvtrksum,	
		id_etaHist,
		id_etaCutHist,
		id_etaHGNHist,
		id_sumtksum_rlo,
		id_sumtksum_rhi,
		id_eRnum,
		id_eRnum_b2p5,
		id_eRnum_nowt,
		id_eRnum_b2p5_nowt,
		id_ePtnum, 
		id_etksumnum,
		id_eRnumf,
                id_eRnumf_b2p5,
                id_ePtnumf,
                id_etksumnumf,
		id_eRden,
		id_eRden_b2p5,
		id_eRden_nowt,
		id_eRden_b2p5_nowt,
		id_ePtden,  
		id_etksumden,
		id_debug1,
		id_debug2,
		id_radflux,
		id_radflux_nowt,
		id_radfluxcoarse,
		id_fluxcomp,
		id_matchdR,
		id_gflux,
		id_gflux_eta,
		id_gflux_nowt,
		id_pc_chi2ndof,
    	    numTH1Hist};
       
	enum th2d_ids{
		id_xyHist,
		id_xywideHist,
		id_rphiHist,
		id_ndof_pcReta,
		id_ndof_pcRpt,
		id_ndof_pcHeta,
                id_ndof_pcHpt,
		id_effptr_num,
		id_effptr_den,
		id_effptr_num_nowt,
		id_effptr_den_nowt,
		id_pt_leadfound,
		id_pt_subfound,
		id_pt_leadlost,
		id_pt_sublost,
		id_pt_leadqual,
		id_pt_subqual,
		id_pt_shared,
		id_pt_rerr,
		id_r_rerr,
                id_reta_pc,
                id_reta_effN,
		id_reta_effD,
		id_reta_effD_monobin,		
		id_reta_effN_nowt,
        id_reta_effD_nowt,
		id_reta_ng,
		id_r_phi_p6,
		id_effrphiN,
		id_effrphiD,
		id_effrphiN_nowt,
		id_effrphiD_nowt,
		id_nchi2_r,
		id_nchi2_dr,
		id_gfluxE,
		id_effER_den,	
	    numTH2Hist};
