package org.hps.recon.utils;

import hep.aida.IAnalysisFactory;
import hep.aida.IHistogram1D;
import hep.aida.IHistogram2D;
import hep.aida.IHistogramFactory;
import hep.aida.ITree;
import hep.aida.ref.rootwriter.RootFileStore;

import hep.physics.vec.Hep3Vector;
import hep.physics.vec.BasicHep3Vector;
//import hep.physics.matrix.SymmetricMatrix;

import java.util.List;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
//import java.util.HashSet;
import java.io.IOException;
import java.util.Set;


import org.lcsim.event.LCRelation;
import org.lcsim.event.RelationalTable;
//import org.lcsim.event.base.BaseLCRelation;
import org.lcsim.event.base.BaseRelationalTable;

import org.lcsim.event.MCParticle;
import org.lcsim.event.Track;
//import org.lcsim.event.TrackState;
//import org.lcsim.event.base.BaseTrack;
//import org.lcsim.event.base.BaseTrackState;
import org.hps.recon.tracking.TrackUtils;
//import org.hps.recon.utils.TrackTruthMatching_new;
//import org.hps.recon.utils.TrackTruthMatching;
//import org.lcsim.fit.helicaltrack.HelicalTrackFit;

import org.lcsim.event.SimTrackerHit;
import org.lcsim.util.Driver;
import org.lcsim.geometry.Detector;
import org.lcsim.event.EventHeader;
import org.lcsim.event.RawTrackerHit;

//import org.lcsim.lcio.LCIOConstants;

/**
 * This driver creates an MCParticle relation to be persisted for each track collection
 * It also saves a TruthTrack
 */
public class TrackToMCParticleRelationsDriver_alic extends Driver {
    
    //Collection Names
    private String trackCollectionName = "GBLTracks";
    
    //If the tracks are kalman tracks
    private boolean kalmanTracks     = true;

    private double bfield;
    private double bfield_y;

    private boolean debug = false;
    private boolean saveTruthTracks = true;

    //Single MCP must leave n hits on Track to be matched to Track
    private int nGoodHitsRequired = 8;
    //Number of hits required on Track
    private int nHitsRequired = 10;

    //(bestMCP_Nhits/Total_Nhits) must be > purity cut    
    //Where "hit" is counted only once per layer
    private double purityCut = 0.8;

    //Define Plot tools
    private ITree tree;
    private IHistogramFactory histogramFactory;
    private Map<String, IHistogram1D> plots1D;
    private Map<String, IHistogram2D> plots2D;
    boolean enablePlots = false;    

    public void setEnablePlots(boolean input){
        this.enablePlots = input;
        if(enablePlots)
            this.bookHistograms();
    }

    public void saveHistograms() {
        String rootFile = String.format("TrackToMCParticleRelations.root");
        RootFileStore store = new RootFileStore(rootFile);
        try {
            store.open();
            store.add(tree);
            store.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    //Book Histograms
    public void bookHistograms(){
        plots1D = new HashMap<String, IHistogram1D>();
        plots2D = new HashMap<String, IHistogram2D>();
        tree = IAnalysisFactory.create().createTreeFactory().create();
        histogramFactory = IAnalysisFactory.create().createHistogramFactory(tree);

//Plots for all Truth Matched Tracks, regradless of purity

        //ele
        plots1D.put(String.format("ele_reco_track_momentum"),
                histogramFactory.createHistogram1D(String.format("ele_reco_track_momentum"), 800, 0, 8));
        //pos
        plots1D.put(String.format("pos_reco_track_momentum"),
                histogramFactory.createHistogram1D(String.format("pos_reco_track_momentum"), 800, 0, 8));

        //ele
        plots1D.put(String.format("ele_n_mcps_on_track"), histogramFactory.createHistogram1D(String.format("ele_n_mcps_on_track"), 30, 0, 30));

        plots1D.put(String.format("ele_n_hits"), histogramFactory.createHistogram1D(String.format("ele_n_hits"), 15, 0, 15));

        plots1D.put(String.format("ele_n_goodhits"), histogramFactory.createHistogram1D(String.format("ele_n_goodhits"), 15, 0, 15));

        plots1D.put(String.format("ele_n_badhits"), histogramFactory.createHistogram1D(String.format("ele_n_badhits"), 15, 0, 15));

        plots1D.put(String.format("ele_purity"), histogramFactory.createHistogram1D(String.format("ele_purity"), 500, 0, 1));

        plots1D.put(String.format("ele_layers_hit"), histogramFactory.createHistogram1D(String.format("ele_layers_hit"), 15, 0, 15));

        plots1D.put(String.format("ele_n_mcps_on_layer"), histogramFactory.createHistogram1D(String.format("ele_n_mcps_on_layer"), 20, 0, 20));

        plots1D.put(String.format("ele_n_striphits_on_layer"), 
                histogramFactory.createHistogram1D(String.format("ele_n_striphits_on_layer"), 20, 0, 20));

        plots1D.put(String.format("ele_n_mcps_on_striphit"), 
                histogramFactory.createHistogram1D(String.format("ele_n_mcps_on_striphit"), 10, 0, 10));

        plots2D.put(String.format("ele_n_mcps_per_layer"), 
                histogramFactory.createHistogram2D(String.format("ele_n_mcps_per_layer"), 20, 0, 20, 15, 0, 15));

        plots2D.put(String.format("ele_n_striphits_per_layer"),
                histogramFactory.createHistogram2D(String.format("ele_n_striphits_per_layer"), 20, 0, 20, 15, 0, 15));

        plots2D.put(String.format("ele_n_mcps_on_layer_striphits"),
                histogramFactory.createHistogram2D(String.format("ele_n_mcps_on_layer_striphits"), 20, 0, 20, 15, 0, 15));

        plots2D.put(String.format("ele_truth_track_v_matched_mcp_momentum"),
                histogramFactory.createHistogram2D(String.format("ele_truth_track_v_matched_mcp_momentum"), 1600, -8, 8, 1600, -8, 8));

        plots1D.put(String.format("ele_truth_track_momentum"),
                histogramFactory.createHistogram1D(String.format("ele_truth_track_momentum"), 800, 0, 8));

        plots1D.put(String.format("ele_truth_track_mcp_momentum"),
                histogramFactory.createHistogram1D(String.format("ele_truth_track_mcp_momentum"), 800, 0, 8));

        plots1D.put(String.format("ele_truth_track_tanlambda"),
                histogramFactory.createHistogram1D(String.format("ele_truth_track_tanlambda"), 1000, -1, 1));

        plots2D.put(String.format("ele_truth_tracks_trackP_v_nHits_on_track"),
                histogramFactory.createHistogram2D(String.format("ele_truth_tracks_trackP_v_nHits_on_track"), 800, 0, 800, 15, 0, 15));

        plots1D.put(String.format("ele_truth_duplicate_tracks_tanlambda"),
                histogramFactory.createHistogram1D(String.format("ele_truth_duplicate_tracks_tanlambda"), 1000, -1, 1));

        plots1D.put(String.format("ele_truth_duplicate_tracks_momentum"),
                histogramFactory.createHistogram1D(String.format("ele_truth_duplicate_tracks_momentum"), 800, 0, 8));
        //pos
        plots1D.put(String.format("pos_n_mcps_on_track"), histogramFactory.createHistogram1D(String.format("pos_n_mcps_on_track"), 30, 0, 30));

        plots1D.put(String.format("pos_n_hits"), histogramFactory.createHistogram1D(String.format("pos_n_hits"), 15, 0, 15));

        plots1D.put(String.format("pos_n_goodhits"), histogramFactory.createHistogram1D(String.format("pos_n_goodhits"), 15, 0, 15));

        plots1D.put(String.format("pos_n_badhits"), histogramFactory.createHistogram1D(String.format("pos_n_badhits"), 15, 0, 15));

        plots1D.put(String.format("pos_purity"), histogramFactory.createHistogram1D(String.format("pos_purity"), 500, 0, 1));

        plots1D.put(String.format("pos_layers_hit"), histogramFactory.createHistogram1D(String.format("pos_layers_hit"), 15, 0, 15));

        plots1D.put(String.format("pos_n_mcps_on_layer"), histogramFactory.createHistogram1D(String.format("pos_n_mcps_on_layer"), 20, 0, 20));

        plots1D.put(String.format("pos_n_striphits_on_layer"), 
                histogramFactory.createHistogram1D(String.format("pos_n_striphits_on_layer"), 20, 0, 20));

        plots1D.put(String.format("pos_n_mcps_on_striphit"), 
                histogramFactory.createHistogram1D(String.format("pos_n_mcps_on_striphit"), 10, 0, 10));

        plots2D.put(String.format("pos_n_mcps_per_layer"), 
                histogramFactory.createHistogram2D(String.format("pos_n_mcps_per_layer"), 20, 0, 20, 15, 0, 15));

        plots2D.put(String.format("pos_n_striphits_per_layer"),
                histogramFactory.createHistogram2D(String.format("pos_n_striphits_per_layer"), 20, 0, 20, 15, 0, 15));

        plots2D.put(String.format("pos_n_mcps_on_layer_striphits"),
                histogramFactory.createHistogram2D(String.format("pos_n_mcps_on_layer_striphits"), 20, 0, 20, 15, 0, 15));

        plots2D.put(String.format("pos_truth_track_v_matched_mcp_momentum"),
                histogramFactory.createHistogram2D(String.format("pos_truth_track_v_matched_mcp_momentum"), 1600, -8, 8, 1600, -8, 8));

        plots1D.put(String.format("pos_truth_track_momentum"),
                histogramFactory.createHistogram1D(String.format("pos_truth_track_momentum"), 800, 0, 8));

        plots1D.put(String.format("pos_truth_track_mcp_momentum"),
                histogramFactory.createHistogram1D(String.format("pos_truth_track_mcp_momentum"), 800, 0, 8));

        plots1D.put(String.format("pos_truth_track_tanlambda"),
                histogramFactory.createHistogram1D(String.format("pos_truth_track_tanlambda"), 1000, -1, 1));

        plots2D.put(String.format("pos_truth_tracks_trackP_v_nHits_on_track"),
                histogramFactory.createHistogram2D(String.format("pos_truth_tracks_trackP_v_nHits_on_track"), 800, 0, 800, 15, 0, 15));

        plots1D.put(String.format("pos_truth_duplicate_tracks_tanlambda"),
                histogramFactory.createHistogram1D(String.format("pos_truth_duplicate_tracks_tanlambda"), 1000, -1, 1));

        plots1D.put(String.format("pos_truth_duplicate_tracks_momentum"),
                histogramFactory.createHistogram1D(String.format("pos_truth_duplicate_tracks_momentum"), 800, 0, 8));

//Plots for "Real" Tracks that pass purity cut

        plots1D.put(String.format("ele_real_track_n_mcps_on_track"),
                histogramFactory.createHistogram1D(String.format("ele_real_track_n_mcps_on_track"), 30, 0, 30));

        plots1D.put(String.format("ele_real_track_n_hits"),
                histogramFactory.createHistogram1D(String.format("ele_real_track_n_hits"), 15, 0, 15));

        plots1D.put(String.format("ele_real_track_n_goodhits"),
                histogramFactory.createHistogram1D(String.format("ele_real_track_n_goodhits"), 15, 0, 15));

        plots1D.put(String.format("ele_real_track_n_badhits"),
                histogramFactory.createHistogram1D(String.format("ele_real_track_n_badhits"), 15, 0, 15));

        plots1D.put(String.format("ele_real_track_purity"),
                histogramFactory.createHistogram1D(String.format("ele_real_track_purity"), 500, 0, 1));

        plots1D.put(String.format("ele_real_track_layers_hit"),
                histogramFactory.createHistogram1D(String.format("ele_real_track_layers_hit"), 15, 0, 15));

        plots1D.put(String.format("ele_real_track_n_mcps_on_layer"),
                histogramFactory.createHistogram1D(String.format("ele_real_track_n_mcps_on_layer"), 20, 0, 20));

        plots1D.put(String.format("ele_real_track_n_striphits_on_layer"), 
                histogramFactory.createHistogram1D(String.format("ele_real_track_n_striphits_on_layer"), 20, 0, 20));

        plots1D.put(String.format("ele_real_track_n_mcps_on_striphit"), 
                histogramFactory.createHistogram1D(String.format("ele_real_track_n_mcps_on_striphit"), 10, 0, 10));

        plots2D.put(String.format("ele_real_track_n_mcps_per_layer"), 
                histogramFactory.createHistogram2D(String.format("ele_real_track_n_mcps_per_layer"), 20, 0, 20, 15, 0, 15));

        plots2D.put(String.format("ele_real_track_n_striphits_per_layer"),
                histogramFactory.createHistogram2D(String.format("ele_real_track_n_striphits_per_layer"), 20, 0, 20, 15, 0, 15));

        plots2D.put(String.format("ele_real_track_n_mcps_on_layer_striphits"),
                histogramFactory.createHistogram2D(String.format("ele_real_track_n_mcps_on_layer_striphits"), 20, 0, 20, 15, 0, 15));

        plots2D.put(String.format("ele_real_track_v_matched_mcp_momentum"),
                histogramFactory.createHistogram2D(String.format("ele_real_track_v_matched_mcp_momentum"), 1600, -8, 8, 1600, -8, 8));

        plots1D.put(String.format("ele_real_track_momentum"),
                histogramFactory.createHistogram1D(String.format("ele_real_track_momentum"), 800, 0, 8));

        plots1D.put(String.format("ele_real_track_mcp_momentum"),
                histogramFactory.createHistogram1D(String.format("ele_real_track_mcp_momentum"), 800, 0, 8));

        plots1D.put(String.format("ele_real_track_tanlambda"),
                histogramFactory.createHistogram1D(String.format("ele_real_track_tanlambda"), 1000, -1, 1));

        plots1D.put(String.format("ele_real_duplicate_tracks_tanlambda"),
                histogramFactory.createHistogram1D(String.format("ele_real_duplicate_tracks_tanlambda"), 1000, -1, 1));

        plots1D.put(String.format("ele_real_duplicate_tracks_momentum"),
                histogramFactory.createHistogram1D(String.format("ele_real_duplicate_tracks_momentum"), 800, 0, 8));

        plots1D.put(String.format("ele_real_duplicate_tracks_mcp_momentum"),
                histogramFactory.createHistogram1D(String.format("ele_real_duplicate_tracks_mcp_momentum"), 800, 0, 8));
        //pos

        plots1D.put(String.format("pos_real_track_n_mcps_on_track"),
                histogramFactory.createHistogram1D(String.format("pos_real_track_n_mcps_on_track"), 30, 0, 30));

        plots1D.put(String.format("pos_real_track_n_hits"),
                histogramFactory.createHistogram1D(String.format("pos_real_track_n_hits"), 15, 0, 15));

        plots1D.put(String.format("pos_real_track_n_goodhits"),
                histogramFactory.createHistogram1D(String.format("pos_real_track_n_goodhits"), 15, 0, 15));

        plots1D.put(String.format("pos_real_track_n_badhits"),
                histogramFactory.createHistogram1D(String.format("pos_real_track_n_badhits"), 15, 0, 15));

        plots1D.put(String.format("pos_real_track_purity"),
                histogramFactory.createHistogram1D(String.format("pos_real_track_purity"), 500, 0, 1));

        plots1D.put(String.format("pos_real_track_layers_hit"),
                histogramFactory.createHistogram1D(String.format("pos_real_track_layers_hit"), 15, 0, 15));

        plots1D.put(String.format("pos_real_track_n_mcps_on_layer"),
                histogramFactory.createHistogram1D(String.format("pos_real_track_n_mcps_on_layer"), 20, 0, 20));

        plots1D.put(String.format("pos_real_track_n_striphits_on_layer"), 
                histogramFactory.createHistogram1D(String.format("pos_real_track_n_striphits_on_layer"), 20, 0, 20));

        plots1D.put(String.format("pos_real_track_n_mcps_on_striphit"), 
                histogramFactory.createHistogram1D(String.format("pos_real_track_n_mcps_on_striphit"), 10, 0, 10));

        plots2D.put(String.format("pos_real_track_n_mcps_per_layer"), 
                histogramFactory.createHistogram2D(String.format("pos_real_track_n_mcps_per_layer"), 20, 0, 20, 15, 0, 15));

        plots2D.put(String.format("pos_real_track_n_striphits_per_layer"),
                histogramFactory.createHistogram2D(String.format("pos_real_track_n_striphits_per_layer"), 20, 0, 20, 15, 0, 15));

        plots2D.put(String.format("pos_real_track_n_mcps_on_layer_striphits"),
                histogramFactory.createHistogram2D(String.format("pos_real_track_n_mcps_on_layer_striphits"), 20, 0, 20, 15, 0, 15));

        plots2D.put(String.format("pos_real_track_v_matched_mcp_momentum"),
                histogramFactory.createHistogram2D(String.format("pos_real_track_v_matched_mcp_momentum"), 1600, -8, 8, 1600, -8, 8));

        plots1D.put(String.format("pos_real_track_momentum"),
                histogramFactory.createHistogram1D(String.format("pos_real_track_momentum"), 800, 0, 8));

        plots1D.put(String.format("pos_real_track_mcp_momentum"),
                histogramFactory.createHistogram1D(String.format("pos_real_track_mcp_momentum"), 800, 0, 8));

        plots1D.put(String.format("pos_real_track_tanlambda"),
                histogramFactory.createHistogram1D(String.format("pos_real_track_tanlambda"), 1000, -1, 1));

        plots1D.put(String.format("pos_real_duplicate_tracks_tanlambda"),
                histogramFactory.createHistogram1D(String.format("pos_real_duplicate_tracks_tanlambda"), 1000, -1, 1));

        plots1D.put(String.format("pos_real_duplicate_tracks_momentum"),
                histogramFactory.createHistogram1D(String.format("pos_real_duplicate_tracks_momentum"), 800, 0, 8));

        plots1D.put(String.format("pos_real_duplicate_tracks_mcp_momentum"),
                histogramFactory.createHistogram1D(String.format("pos_real_duplicate_tracks_mcp_momentum"), 800, 0, 8));

//Plots for "Fake" Tracks that fail purity cut

        //ele
        plots1D.put(String.format("ele_fake_track_n_mcps_on_track"),
                histogramFactory.createHistogram1D(String.format("ele_fake_track_n_mcps_on_track"), 30, 0, 30));

        plots1D.put(String.format("ele_fake_track_n_hits"),
                histogramFactory.createHistogram1D(String.format("ele_fake_track_n_hits"), 15, 0, 15));

        plots1D.put(String.format("ele_fake_track_n_goodhits"),
                histogramFactory.createHistogram1D(String.format("ele_fake_track_n_goodhits"), 15, 0, 15));

        plots1D.put(String.format("ele_fake_track_n_badhits"),
                histogramFactory.createHistogram1D(String.format("ele_fake_track_n_badhits"), 15, 0, 15));

        plots1D.put(String.format("ele_fake_track_purity"),
                histogramFactory.createHistogram1D(String.format("ele_fake_track_purity"), 500, 0, 1));

        plots1D.put(String.format("ele_fake_track_layers_hit"),
                histogramFactory.createHistogram1D(String.format("ele_fake_track_layers_hit"), 15, 0, 15));

        plots1D.put(String.format("ele_fake_track_n_mcps_on_layer"),
                histogramFactory.createHistogram1D(String.format("ele_fake_track_n_mcps_on_layer"), 20, 0, 20));

        plots1D.put(String.format("ele_fake_track_n_striphits_on_layer"), 
                histogramFactory.createHistogram1D(String.format("ele_fake_track_n_striphits_on_layer"), 20, 0, 20));

        plots1D.put(String.format("ele_fake_track_n_mcps_on_striphit"), 
                histogramFactory.createHistogram1D(String.format("ele_fake_track_n_mcps_on_striphit"), 10, 0, 10));

        plots2D.put(String.format("ele_fake_track_n_mcps_per_layer"), 
                histogramFactory.createHistogram2D(String.format("ele_fake_track_n_mcps_per_layer"), 20, 0, 20, 15, 0, 15));

        plots2D.put(String.format("ele_fake_track_n_striphits_per_layer"),
                histogramFactory.createHistogram2D(String.format("ele_fake_track_n_striphits_per_layer"), 20, 0, 20, 15, 0, 15));

        plots2D.put(String.format("ele_fake_track_n_mcps_on_layer_striphits"),
                histogramFactory.createHistogram2D(String.format("ele_fake_track_n_mcps_on_layer_striphits"), 20, 0, 20, 15, 0, 15));

        plots2D.put(String.format("ele_fake_track_v_matched_mcp_momentum"),
                histogramFactory.createHistogram2D(String.format("ele_fake_track_v_matched_mcp_momentum"), 1600, -8, 8, 1600, -8, 8));

        plots1D.put(String.format("ele_fake_track_momentum"),
                histogramFactory.createHistogram1D(String.format("ele_fake_track_momentum"), 800, 0, 8));

        plots1D.put(String.format("ele_fake_track_mcp_momentum"),
                histogramFactory.createHistogram1D(String.format("ele_fake_track_mcp_momentum"), 800, 0, 8));

        plots1D.put(String.format("ele_fake_track_tanlambda"),
                histogramFactory.createHistogram1D(String.format("ele_fake_track_tanlambda"), 1000, -1, 1));

        plots1D.put(String.format("ele_fake_duplicate_tracks_tanlambda"),
                histogramFactory.createHistogram1D(String.format("ele_fake_duplicate_tracks_tanlambda"), 1000, -1, 1));

        plots1D.put(String.format("ele_fake_duplicate_tracks_momentum"),
                histogramFactory.createHistogram1D(String.format("ele_fake_duplicate_tracks_momentum"), 800, 0, 8));

        plots1D.put(String.format("ele_fake_duplicate_tracks_mcp_momentum"),
                histogramFactory.createHistogram1D(String.format("ele_fake_duplicate_tracks_mcp_momentum"), 800, 0, 8));
        //pos
        plots1D.put(String.format("pos_fake_track_n_mcps_on_track"),
                histogramFactory.createHistogram1D(String.format("pos_fake_track_n_mcps_on_track"), 30, 0, 30));

        plots1D.put(String.format("pos_fake_track_n_hits"),
                histogramFactory.createHistogram1D(String.format("pos_fake_track_n_hits"), 15, 0, 15));

        plots1D.put(String.format("pos_fake_track_n_goodhits"),
                histogramFactory.createHistogram1D(String.format("pos_fake_track_n_goodhits"), 15, 0, 15));

        plots1D.put(String.format("pos_fake_track_n_badhits"),
                histogramFactory.createHistogram1D(String.format("pos_fake_track_n_badhits"), 15, 0, 15));

        plots1D.put(String.format("pos_fake_track_purity"),
                histogramFactory.createHistogram1D(String.format("pos_fake_track_purity"), 500, 0, 1));

        plots1D.put(String.format("pos_fake_track_layers_hit"),
                histogramFactory.createHistogram1D(String.format("pos_fake_track_layers_hit"), 15, 0, 15));

        plots1D.put(String.format("pos_fake_track_n_mcps_on_layer"),
                histogramFactory.createHistogram1D(String.format("pos_fake_track_n_mcps_on_layer"), 20, 0, 20));

        plots1D.put(String.format("pos_fake_track_n_striphits_on_layer"), 
                histogramFactory.createHistogram1D(String.format("pos_fake_track_n_striphits_on_layer"), 20, 0, 20));

        plots1D.put(String.format("pos_fake_track_n_mcps_on_striphit"), 
                histogramFactory.createHistogram1D(String.format("pos_fake_track_n_mcps_on_striphit"), 10, 0, 10));

        plots2D.put(String.format("pos_fake_track_n_mcps_per_layer"), 
                histogramFactory.createHistogram2D(String.format("pos_fake_track_n_mcps_per_layer"), 20, 0, 20, 15, 0, 15));

        plots2D.put(String.format("pos_fake_track_n_striphits_per_layer"),
                histogramFactory.createHistogram2D(String.format("pos_fake_track_n_striphits_per_layer"), 20, 0, 20, 15, 0, 15));

        plots2D.put(String.format("pos_fake_track_n_mcps_on_layer_striphits"),
                histogramFactory.createHistogram2D(String.format("pos_fake_track_n_mcps_on_layer_striphits"), 20, 0, 20, 15, 0, 15));

        plots2D.put(String.format("pos_fake_track_v_matched_mcp_momentum"),
                histogramFactory.createHistogram2D(String.format("pos_fake_track_v_matched_mcp_momentum"), 1600, -8, 8, 1600, -8, 8));

        plots1D.put(String.format("pos_fake_track_momentum"),
                histogramFactory.createHistogram1D(String.format("pos_fake_track_momentum"), 800, 0, 8));

        plots1D.put(String.format("pos_fake_track_mcp_momentum"),
                histogramFactory.createHistogram1D(String.format("pos_fake_track_mcp_momentum"), 800, 0, 8));

        plots1D.put(String.format("pos_fake_track_tanlambda"),
                histogramFactory.createHistogram1D(String.format("pos_fake_track_tanlambda"), 1000, -1, 1));

        plots1D.put(String.format("pos_fake_duplicate_tracks_tanlambda"),
                histogramFactory.createHistogram1D(String.format("pos_fake_duplicate_tracks_tanlambda"), 1000, -1, 1));

        plots1D.put(String.format("pos_fake_duplicate_tracks_momentum"),
                histogramFactory.createHistogram1D(String.format("pos_fake_duplicate_tracks_momentum"), 800, 0, 8));

        plots1D.put(String.format("pos_fake_duplicate_tracks_mcp_momentum"),
                histogramFactory.createHistogram1D(String.format("pos_fake_duplicate_tracks_mcp_momentum"), 800, 0, 8));

//
        plots2D.put(String.format("old_truth_track_v_matched_mcp_momentum"),
                histogramFactory.createHistogram2D(String.format("old_truth_track_v_matched_mcp_momentum"), 1600, -8, 8, 1600, -8, 8));

        //Trackable MCP Plots
        plots1D.put(String.format("mcp_px"), histogramFactory.createHistogram1D(String.format("mcp_px"), 1000, -10, 10));

        plots1D.put(String.format("mcp_py"), histogramFactory.createHistogram1D(String.format("mcp_py"), 1000, -10, 10));
        
        plots1D.put(String.format("mcp_pz"), histogramFactory.createHistogram1D(String.format("mcp_pz"), 1000, -10, 10));

        plots2D.put(String.format("mcp_momentum_v_px"),
                histogramFactory.createHistogram2D(String.format("mcp_momentum_v_px"), 800, 0, 8, 1000, -10, 10));
        
        plots2D.put(String.format("mcp_momentum_v_py"),
                histogramFactory.createHistogram2D(String.format("mcp_momentum_v_py"), 800, 0, 8, 1000, -10, 10));
        
        plots2D.put(String.format("mcp_momentum_v_pz"),
                histogramFactory.createHistogram2D(String.format("mcp_momentum_v_pz"), 800, 0, 8, 1000, -10, 10));

        //ele
        plots1D.put(String.format("ele_trackable_mcp_nSimTrackerHits"),
                histogramFactory.createHistogram1D(String.format("ele_trackable_mcp_nSimTrackerHits"), 15, 0, 15));
        
        plots2D.put(String.format("ele_trackable_mcp_momentum_v_nSimTrackerHits"),
                histogramFactory.createHistogram2D(String.format("ele_trackable_mcp_momentum_v_nSimTrackerHits"), 800, 0, 8, 15, 0, 15));

        plots1D.put(String.format("ele_trackable_mcp_momentum"),
                histogramFactory.createHistogram1D(String.format("ele_trackable_mcp_momentum"), 800, 0, 8));
        //pos
        plots1D.put(String.format("pos_trackable_mcp_nSimTrackerHits"),
                histogramFactory.createHistogram1D(String.format("pos_trackable_mcp_nSimTrackerHits"), 15, 0, 15));
        
        plots2D.put(String.format("pos_trackable_mcp_momentum_v_nSimTrackerHits"),
                histogramFactory.createHistogram2D(String.format("pos_trackable_mcp_momentum_v_nSimTrackerHits"), 800, 0, 8, 15, 0, 15));

        plots1D.put(String.format("pos_trackable_mcp_momentum"),
                histogramFactory.createHistogram1D(String.format("pos_trackable_mcp_momentum"), 800, 0, 8));
    }

    public void endOfData() {
        if(enablePlots)
            saveHistograms();
    }

    public void setTrackCollectionName(String trackCollectionName) {
        this.trackCollectionName = trackCollectionName;
    }

    public void setKalmanTracks(boolean val) {
        kalmanTracks = val;
    }

    public void setDebug(boolean val) {
        debug = val;
    }
    
    public void setSaveTruthTracks(boolean val) {
        saveTruthTracks = val;
    }
    
    /*
      public void setTrackCollectionNames(String [] trackCollectionNames) {
      this.trackCollectionNames = trackCollectionNames;
      }
    */
    
    @Override
    protected void detectorChanged(Detector detector) {
        Hep3Vector fieldInTracker = TrackUtils.getBField(detector);
        bfield = Math.abs(fieldInTracker.y());
        bfield_y = fieldInTracker.y();
    }
    
    @Override 
    protected void process(EventHeader event) {
        
        //Retrieve track collection
        List<Track> trackCollection = new ArrayList<Track>();
        if (trackCollectionName != null) {
            if (event.hasCollection(Track.class, trackCollectionName)) {
                trackCollection = event.get(Track.class, trackCollectionName);
            }
        }
        
        //Retrieve rawhits to mc
        
        RelationalTable rawtomc = new BaseRelationalTable(RelationalTable.Mode.MANY_TO_MANY, RelationalTable.Weighting.UNWEIGHTED);
        if (event.hasCollection(LCRelation.class, "SVTTrueHitRelations")) {
            List<LCRelation> trueHitRelations = event.get(LCRelation.class, "SVTTrueHitRelations");
            for (LCRelation relation : trueHitRelations)
                if (relation != null && relation.getFrom() != null && relation.getTo() != null)
                    rawtomc.add(relation.getFrom(), relation.getTo());
        }
        else
            return;

        //Retrieve all simulated hits
        String MCHitInputCollectionName = "TrackerHits";
        List<SimTrackerHit> allsimhits = event.get(SimTrackerHit.class, MCHitInputCollectionName);

        //Retrieve all MCPs
        List<MCParticle> allmcps = event.get(MCParticle.class, "MCParticle");
        //Check for Trackable MCPs
        Map<MCParticle, Map<Integer, List<SimTrackerHit>>> trackableMCPMap = getTrackableMCPs(allmcps, allsimhits, rawtomc,this.nGoodHitsRequired);


        //MCParticleRelations
        List<LCRelation> trackToMCParticleRelations    =  new ArrayList<LCRelation>();
        
        //Truth Tracks and Relations
        List<LCRelation> trackToTruthTrackRelations    =  new ArrayList<LCRelation>();
        List<Track>      truthTrackCollection          =  new ArrayList<Track>();
        
        List<MCParticle> mcps = new ArrayList<MCParticle>();
        for (Track track : trackCollection) {
            boolean realTrack = true;
            int charge = -1* (int)Math.signum(track.getTrackStates().get(0).getOmega());
            double trackPmag = new BasicHep3Vector(track.getTrackStates().get(0).getMomentum()).magnitude();
            double tanlambda = track.getTrackStates().get(0).getTanLambda();

            //Check match quality of old Track->MCP matching
            TrackTruthMatching ttm = new TrackTruthMatching(track, rawtomc, allsimhits, false);
            MCParticle mcp_matt = null;
            if(ttm != null)
                mcp_matt = ttm.getMCParticle();
            if(mcp_matt != null){
                plots2D.get("old_truth_track_v_matched_mcp_momentum").fill(charge*trackPmag,((int)mcp_matt.getCharge())*mcp_matt.getMomentum().magnitude());
            }

            //Use New TrackTruthMatching tool to match Track -> MCP
            TrackTruthMatching_new tt = new TrackTruthMatching_new(track, rawtomc, 0.0, 0);

            //Check nHits on Track
            //If nHits < required, dont analyze this Track
            double nHits = tt.getNHits();
            if(nHits < this.nHitsRequired)
                continue;

            //Fill plot with all Reco Tracks momentum
            if(charge < 0){
                plots1D.get("ele_reco_track_momentum").fill(trackPmag);
            }
            else{
                plots1D.get("pos_reco_track_momentum").fill(trackPmag);
            }

            //get Track MCP if it has one
            MCParticle mcp = tt.getMCParticle();

            //If Track is matched to a MCP, continue analysis
            if(mcp != null){
                //Get purity and nGoodHits to check Track quality
                double purity = tt.getPurity();
                double nGoodHits = tt.getNGoodHits();

                //First fill plots for all cases where a Track is matched to a
                //MCP, regardless of match purity.
                if(charge < 0){
                    plots2D.get("ele_truth_tracks_trackP_v_nHits_on_track").fill(trackPmag, nHits);
                    plots1D.get("ele_truth_track_momentum").fill(trackPmag);
                    plots1D.get("ele_truth_track_tanlambda").fill(tanlambda);
                    if(mcps.contains(mcp)){
                        plots1D.get("ele_truth_duplicate_tracks_momentum").fill(trackPmag);
                        plots1D.get("ele_truth_duplicate_tracks_tanlambda").fill(tanlambda);
                    }
                    plots2D.get("ele_truth_track_v_matched_mcp_momentum").fill(charge*trackPmag,((int)mcp.getCharge())*mcp.getMomentum().magnitude());
                    plots1D.get("ele_n_mcps_on_track").fill(tt.getLayerHitsForAllMCPs().size());
                    plots1D.get("ele_n_hits").fill(tt.getNHits());
                    plots1D.get("ele_purity").fill(tt.getPurity());
                    plots1D.get("ele_n_goodhits").fill(tt.getNGoodHits());
                    plots1D.get("ele_n_badhits").fill(tt.getNBadHits());
                    for(Integer layer : tt.getLayersOnTrack()){
                        plots1D.get("ele_layers_hit").fill(layer);
                        plots1D.get("ele_n_mcps_on_layer").fill(tt.getMCPsOnLayer(layer).size());
                        plots2D.get("ele_n_mcps_per_layer").fill(layer,tt.getMCPsOnLayer(layer).size());
                        plots1D.get("ele_n_striphits_on_layer").fill(tt.getStripHitsOnLayer(layer).size());
                        plots2D.get("ele_n_striphits_per_layer").fill(layer,tt.getStripHitsOnLayer(layer).size());
                        for(RawTrackerHit rawhit : tt.getStripHitsOnLayer(layer)){
                            plots1D.get("ele_n_mcps_on_striphit").fill(tt.getMCPsOnRawTrackerHit(rawhit).size());
                            plots2D.get("ele_n_mcps_on_layer_striphits").fill(layer,tt.getMCPsOnRawTrackerHit(rawhit).size());
                        }
                    }
                }

                else{
                    plots2D.get("pos_truth_tracks_trackP_v_nHits_on_track").fill(trackPmag, nHits);
                    plots1D.get("pos_truth_track_momentum").fill(trackPmag);
                    plots1D.get("pos_truth_track_tanlambda").fill(tanlambda);
                    if(mcps.contains(mcp)){
                        plots1D.get("pos_truth_duplicate_tracks_momentum").fill(trackPmag);
                        plots1D.get("pos_truth_duplicate_tracks_tanlambda").fill(tanlambda);
                    }
                    plots2D.get("pos_truth_track_v_matched_mcp_momentum").fill(charge*trackPmag,((int)mcp.getCharge())*mcp.getMomentum().magnitude());
                    plots1D.get("pos_n_mcps_on_track").fill(tt.getLayerHitsForAllMCPs().size());
                    plots1D.get("pos_n_hits").fill(tt.getNHits());
                    plots1D.get("pos_purity").fill(tt.getPurity());
                    plots1D.get("pos_n_goodhits").fill(tt.getNGoodHits());
                    plots1D.get("pos_n_badhits").fill(tt.getNBadHits());
                    for(Integer layer : tt.getLayersOnTrack()){
                        plots1D.get("pos_layers_hit").fill(layer);
                        plots1D.get("pos_n_mcps_on_layer").fill(tt.getMCPsOnLayer(layer).size());
                        plots2D.get("pos_n_mcps_per_layer").fill(layer,tt.getMCPsOnLayer(layer).size());
                        plots1D.get("pos_n_striphits_on_layer").fill(tt.getStripHitsOnLayer(layer).size());
                        plots2D.get("pos_n_striphits_per_layer").fill(layer,tt.getStripHitsOnLayer(layer).size());
                        for(RawTrackerHit rawhit : tt.getStripHitsOnLayer(layer)){
                            plots1D.get("pos_n_mcps_on_striphit").fill(tt.getMCPsOnRawTrackerHit(rawhit).size());
                            plots2D.get("pos_n_mcps_on_layer_striphits").fill(layer,tt.getMCPsOnRawTrackerHit(rawhit).size());
                        }
                    }
                }


                //If Track passes purity cut, identify Track as real
                if(purity < this.purityCut || nHits < this.nHitsRequired)
                    realTrack = false;

                //If Track fails purity cut, identify Track as fake
                if(!realTrack){

                    if(charge < 0){
                        plots1D.get("ele_fake_track_momentum").fill(trackPmag);
                        plots1D.get("ele_fake_track_tanlambda").fill(tanlambda);
                        if(mcps.contains(mcp)){
                            plots1D.get("ele_fake_duplicate_tracks_momentum").fill(trackPmag);
                            plots1D.get("ele_fake_duplicate_tracks_tanlambda").fill(tanlambda);
                        }
                        plots2D.get("ele_fake_track_v_matched_mcp_momentum").fill(charge*trackPmag,((int)mcp.getCharge())*mcp.getMomentum().magnitude());
                        plots1D.get("ele_fake_track_n_mcps_on_track").fill(tt.getLayerHitsForAllMCPs().size());
                        plots1D.get("ele_fake_track_n_hits").fill(tt.getNHits());
                        plots1D.get("ele_fake_track_purity").fill(tt.getPurity());
                        plots1D.get("ele_fake_track_n_goodhits").fill(tt.getNGoodHits());
                        plots1D.get("ele_fake_track_n_badhits").fill(tt.getNBadHits());
                        for(Integer layer : tt.getLayersOnTrack()){
                            plots1D.get("ele_fake_track_layers_hit").fill(layer);
                            plots1D.get("ele_fake_track_n_mcps_on_layer").fill(tt.getMCPsOnLayer(layer).size());
                            plots2D.get("ele_fake_track_n_mcps_per_layer").fill(layer,tt.getMCPsOnLayer(layer).size());
                            plots1D.get("ele_fake_track_n_striphits_on_layer").fill(tt.getStripHitsOnLayer(layer).size());
                            plots2D.get("ele_fake_track_n_striphits_per_layer").fill(layer,tt.getStripHitsOnLayer(layer).size());
                            for(RawTrackerHit rawhit : tt.getStripHitsOnLayer(layer)){
                                plots1D.get("ele_fake_track_n_mcps_on_striphit").fill(tt.getMCPsOnRawTrackerHit(rawhit).size());
                                plots2D.get("ele_fake_track_n_mcps_on_layer_striphits").fill(layer,tt.getMCPsOnRawTrackerHit(rawhit).size());
                            }
                        }
                    }

                    else{
                        plots1D.get("pos_fake_track_momentum").fill(trackPmag);
                        plots1D.get("pos_fake_track_tanlambda").fill(tanlambda);
                        if(mcps.contains(mcp)){
                            plots1D.get("pos_fake_duplicate_tracks_momentum").fill(trackPmag);
                            plots1D.get("pos_fake_duplicate_tracks_tanlambda").fill(tanlambda);
                        }
                        plots2D.get("pos_fake_track_v_matched_mcp_momentum").fill(charge*trackPmag,((int)mcp.getCharge())*mcp.getMomentum().magnitude());
                        plots1D.get("pos_fake_track_n_mcps_on_track").fill(tt.getLayerHitsForAllMCPs().size());
                        plots1D.get("pos_fake_track_n_hits").fill(tt.getNHits());
                        plots1D.get("pos_fake_track_purity").fill(tt.getPurity());
                        plots1D.get("pos_fake_track_n_goodhits").fill(tt.getNGoodHits());
                        plots1D.get("pos_fake_track_n_badhits").fill(tt.getNBadHits());
                        for(Integer layer : tt.getLayersOnTrack()){
                            plots1D.get("pos_fake_track_layers_hit").fill(layer);
                            plots1D.get("pos_fake_track_n_mcps_on_layer").fill(tt.getMCPsOnLayer(layer).size());
                            plots2D.get("pos_fake_track_n_mcps_per_layer").fill(layer,tt.getMCPsOnLayer(layer).size());
                            plots1D.get("pos_fake_track_n_striphits_on_layer").fill(tt.getStripHitsOnLayer(layer).size());
                            plots2D.get("pos_fake_track_n_striphits_per_layer").fill(layer,tt.getStripHitsOnLayer(layer).size());
                            for(RawTrackerHit rawhit : tt.getStripHitsOnLayer(layer)){
                                plots1D.get("pos_fake_track_n_mcps_on_striphit").fill(tt.getMCPsOnRawTrackerHit(rawhit).size());
                                plots2D.get("pos_fake_track_n_mcps_on_layer_striphits").fill(layer,tt.getMCPsOnRawTrackerHit(rawhit).size());
                            }
                        }
                    }

                    //MCP based plots
                    if(mcp.getCharge() < 0){
                        plots1D.get("ele_fake_track_mcp_momentum").fill(mcp.getMomentum().magnitude());
                        if(mcps.contains(mcp)){
                            plots1D.get("ele_fake_duplicate_tracks_mcp_momentum").fill(mcp.getMomentum().magnitude());
                        }
                    }
                    else{
                        plots1D.get("pos_fake_track_mcp_momentum").fill(mcp.getMomentum().magnitude());
                        if(mcps.contains(mcp)){
                            plots1D.get("pos_fake_duplicate_tracks_mcp_momentum").fill(mcp.getMomentum().magnitude());
                        }

                    }
                }

                //If Track is real
                if(realTrack){

                    if(charge < 0){
                        plots1D.get("ele_real_track_momentum").fill(trackPmag);
                        plots1D.get("ele_real_track_tanlambda").fill(tanlambda);
                        if(mcps.contains(mcp)){
                            plots1D.get("ele_real_duplicate_tracks_momentum").fill(trackPmag);
                            plots1D.get("ele_real_duplicate_tracks_tanlambda").fill(tanlambda);
                        }
                        plots2D.get("ele_real_track_v_matched_mcp_momentum").fill(charge*trackPmag,((int)mcp.getCharge())*mcp.getMomentum().magnitude());
                        plots1D.get("ele_real_track_n_mcps_on_track").fill(tt.getLayerHitsForAllMCPs().size());
                        plots1D.get("ele_real_track_n_hits").fill(tt.getNHits());
                        plots1D.get("ele_real_track_purity").fill(tt.getPurity());
                        plots1D.get("ele_real_track_n_goodhits").fill(tt.getNGoodHits());
                        plots1D.get("ele_real_track_n_badhits").fill(tt.getNBadHits());
                        for(Integer layer : tt.getLayersOnTrack()){
                            plots1D.get("ele_real_track_layers_hit").fill(layer);
                            plots1D.get("ele_real_track_n_mcps_on_layer").fill(tt.getMCPsOnLayer(layer).size());
                            plots2D.get("ele_real_track_n_mcps_per_layer").fill(layer,tt.getMCPsOnLayer(layer).size());
                            plots1D.get("ele_real_track_n_striphits_on_layer").fill(tt.getStripHitsOnLayer(layer).size());
                            plots2D.get("ele_real_track_n_striphits_per_layer").fill(layer,tt.getStripHitsOnLayer(layer).size());
                            for(RawTrackerHit rawhit : tt.getStripHitsOnLayer(layer)){
                                plots1D.get("ele_real_track_n_mcps_on_striphit").fill(tt.getMCPsOnRawTrackerHit(rawhit).size());
                                plots2D.get("ele_real_track_n_mcps_on_layer_striphits").fill(layer,tt.getMCPsOnRawTrackerHit(rawhit).size());
                            }
                        }
                    }


                    else{
                        plots1D.get("pos_real_track_momentum").fill(trackPmag);
                        plots1D.get("pos_real_track_tanlambda").fill(tanlambda);
                        if(mcps.contains(mcp)){
                            plots1D.get("pos_real_duplicate_tracks_momentum").fill(trackPmag);
                            plots1D.get("pos_real_duplicate_tracks_tanlambda").fill(tanlambda);
                        }
                        plots2D.get("pos_real_track_v_matched_mcp_momentum").fill(charge*trackPmag,((int)mcp.getCharge())*mcp.getMomentum().magnitude());
                        plots1D.get("pos_real_track_n_mcps_on_track").fill(tt.getLayerHitsForAllMCPs().size());
                        plots1D.get("pos_real_track_n_hits").fill(tt.getNHits());
                        plots1D.get("pos_real_track_purity").fill(tt.getPurity());
                        plots1D.get("pos_real_track_n_goodhits").fill(tt.getNGoodHits());
                        plots1D.get("pos_real_track_n_badhits").fill(tt.getNBadHits());
                        for(Integer layer : tt.getLayersOnTrack()){
                            plots1D.get("pos_real_track_layers_hit").fill(layer);
                            plots1D.get("pos_real_track_n_mcps_on_layer").fill(tt.getMCPsOnLayer(layer).size());
                            plots2D.get("pos_real_track_n_mcps_per_layer").fill(layer,tt.getMCPsOnLayer(layer).size());
                            plots1D.get("pos_real_track_n_striphits_on_layer").fill(tt.getStripHitsOnLayer(layer).size());
                            plots2D.get("pos_real_track_n_striphits_per_layer").fill(layer,tt.getStripHitsOnLayer(layer).size());
                            for(RawTrackerHit rawhit : tt.getStripHitsOnLayer(layer)){
                                plots1D.get("pos_real_track_n_mcps_on_striphit").fill(tt.getMCPsOnRawTrackerHit(rawhit).size());
                                plots2D.get("pos_real_track_n_mcps_on_layer_striphits").fill(layer,tt.getMCPsOnRawTrackerHit(rawhit).size());
                            }
                        }
                    }

                    //MCP based plots
                    if(mcp.getCharge() < 0){
                        plots1D.get("ele_real_track_mcp_momentum").fill(mcp.getMomentum().magnitude());
                        if(mcps.contains(mcp)){
                            plots1D.get("ele_real_duplicate_tracks_mcp_momentum").fill(mcp.getMomentum().magnitude());
                        }
                    }
                    else{

                        plots1D.get("pos_real_track_mcp_momentum").fill(mcp.getMomentum().magnitude());
                        if(mcps.contains(mcp)){
                            plots1D.get("pos_real_duplicate_tracks_mcp_momentum").fill(mcp.getMomentum().magnitude());
                        }
                    }
                }

                //Add this MCP to a list to check for duplicate Track->MCP
                //matches 
                mcps.add(mcp);
            }

            //No MCP matched to this Track, so Track is fake...
            else{

                if(charge < 0){
                    plots1D.get("ele_fake_track_momentum").fill(trackPmag);
                    plots1D.get("ele_fake_track_tanlambda").fill(tanlambda);
                }
                else{
                    plots1D.get("pos_fake_track_momentum").fill(trackPmag);
                    plots1D.get("pos_fake_track_tanlambda").fill(tanlambda);
                }

                //MCP based plots
                if(mcp.getCharge() < 0){
                    plots1D.get("ele_fake_track_mcp_momentum").fill(mcp.getMomentum().magnitude());
                }
                else{
                    plots1D.get("pos_fake_track_mcp_momentum").fill(mcp.getMomentum().magnitude());
                }

            }

            /*
            //Truth Matching tool
            TrackTruthMatching ttm = new TrackTruthMatching(track, rawtomc, allsimhits, kalmanTracks);
            if (ttm != null) {
                MCParticle mcp = ttm.getMCParticle();
                
                if (mcp != null) {
                    trackToMCParticleRelations.add(new BaseLCRelation(track,mcp));
                    
                    //Hep3Vector origin = new BasicHep3Vector(0.,0.,0.);
                    HelicalTrackFit mcp_htf  = TrackUtils.getHTF(mcp,bfield);
                    HelicalTrackFit trk_htf  = TrackUtils.getHTF(track);
                
                    if (debug) {
                        System.out.println("--------------------");
                        System.out.println(trackCollectionName+" Track:");
                        System.out.printf("d0 %f z0 %f R %f phi %f lambda %s\n", trk_htf.dca(), trk_htf.z0(), trk_htf.R(), trk_htf.phi0(), trk_htf.slope());
                        System.out.printf("Nhits = %d \n",ttm.getNHits());
                        System.out.printf("NGoodHits = %d  purity = %f\n",ttm.getNGoodHits(),ttm.getPurity());
                        
                    }
                    BaseTrack truth_trk  = new BaseTrack();
                    truth_trk.setTrackParameters(mcp_htf.parameters(),bfield);
                    truth_trk.getTrackStates().clear();
                    double[] ref = new double[] { 0., 0., 0. };
                    SymmetricMatrix cov = new SymmetricMatrix(5);
                    TrackState stateIP = new BaseTrackState(mcp_htf.parameters(),ref,cov.asPackedArray(true),TrackState.AtIP,bfield);
                    truth_trk.getTrackStates().add(stateIP);
                    truth_trk.setChisq(-1);
                    truth_trk.setNDF(-1);
                    truth_trk.setFitSuccess(false);
                    truth_trk.setRefPointIsDCA(true);
                    truth_trk.setTrackType(-1);
                    truthTrackCollection.add(truth_trk);
                    trackToTruthTrackRelations.add(new BaseLCRelation(track,truth_trk));
                    
                    
                    if (debug) {
                        double d0    = truth_trk.getTrackStates().get(0).getD0();
                        double z0    = truth_trk.getTrackStates().get(0).getZ0();
                        double C     = truth_trk.getTrackStates().get(0).getOmega();
                        double phi   = truth_trk.getTrackStates().get(0).getPhi();
                        double slope = truth_trk.getTrackStates().get(0).getTanLambda();
                        System.out.printf("TruthTrack \n");
                        System.out.printf("d0 %f z0 %f R %f phi %f lambda %s\n", d0, z0, 1./C, phi, slope);
                        System.out.printf("MCParticle  \n");
                        Hep3Vector pVec = mcp.getMomentum();
                        System.out.printf("ptTrue = %f \n", Math.sqrt(pVec.x()*pVec.x() + pVec.z()*pVec.z()));
                        double mom_param = 2.99792458e-04;
                        double trkMom = trk_htf.R() * bfield * mom_param;
                        System.out.printf("pt = %f \n", trkMom);
                        
                        System.out.println("--------------------");
                    }
                }//mcp not null
                else {
                    System.out.printf("PF::FakeTrack");
                }
            } //ttm not null
            else {
                System.out.printf("Error::TrackTruthMatching returns null \n");
            }
        */
        }
        
        //int flag = 1 << LCIOConstants.TRBIT_HITS;
        //event.put(trackCollectionName+"Truth", truthTrackCollection, Track.class, flag);
        //event.put(trackCollectionName+"ToTruthTrackRelations", trackToTruthTrackRelations, LCRelation.class, 0);
        //event.put(trackCollectionName+"ToMCParticleRelations", trackToMCParticleRelations, LCRelation.class, 0);
    }//closes process

    public Map<MCParticle, Map<Integer, List<SimTrackerHit>>> getTrackableMCPs(List<MCParticle> mcparticles, List<SimTrackerHit> simhits, RelationalTable rawtomc, int NhitsRequired){
        Map<MCParticle, Map<Integer, List<SimTrackerHit>>> trackableMCPs = new HashMap<MCParticle, Map<Integer, List<SimTrackerHit>>>();

        for(MCParticle particle : mcparticles){
            double P = particle.getMomentum().magnitude();
            if(Math.abs(particle.getPDGID()) != 11)
                continue;
            //if(!particle.getParents().contains(622) && !particle.getParents().contains(623) && P < fee) 
                //continue;
            double px = particle.getPX();
            double py = particle.getPY();
            double pz = particle.getPZ();
            plots1D.get("mcp_px").fill(px);
            plots1D.get("mcp_py").fill(py);
            plots1D.get("mcp_pz").fill(pz);
            plots2D.get("mcp_momentum_v_px").fill(P, px);
            plots2D.get("mcp_momentum_v_py").fill(P, py);
            plots2D.get("mcp_momentum_v_pz").fill(P, pz);

            Map<Integer, List<SimTrackerHit>> layerhitsMap = new HashMap<Integer, List<SimTrackerHit>>();
            for(SimTrackerHit simhit : simhits){
                MCParticle simhitmcp = simhit.getMCParticle();
                if(simhitmcp == particle){
                    int layer = simhit.getLayer();
                    Set<RawTrackerHit> rawhits = rawtomc.allTo(simhit);
                    System.out.println("Nrawhits on simhit is " + rawhits.size());
                    if(!layerhitsMap.containsKey(layer)){
                        List<SimTrackerHit> tmp = new ArrayList<SimTrackerHit>();
                        tmp.add(simhit);
                        layerhitsMap.put(layer, tmp); 
                    }
                    else{
                        List<SimTrackerHit> tmp = layerhitsMap.get(layer);
                        tmp.add(simhit);
                        layerhitsMap.put(layer, tmp);
                    }
                }
            }
            if(layerhitsMap.size() < NhitsRequired)
                continue;
            if(particle.getCharge() < 0){
                plots1D.get("ele_trackable_mcp_nSimTrackerHits").fill(layerhitsMap.size());
                plots2D.get("ele_trackable_mcp_momentum_v_nSimTrackerHits").fill(particle.getMomentum().magnitude(), layerhitsMap.size());
                plots1D.get("ele_trackable_mcp_momentum").fill(particle.getMomentum().magnitude());
            }
            else{
                plots1D.get("pos_trackable_mcp_nSimTrackerHits").fill(layerhitsMap.size());
                plots2D.get("pos_trackable_mcp_momentum_v_nSimTrackerHits").fill(particle.getMomentum().magnitude(), layerhitsMap.size());
                plots1D.get("pos_trackable_mcp_momentum").fill(particle.getMomentum().magnitude());
            }

            trackableMCPs.put(particle, layerhitsMap);
        }
        return trackableMCPs;
    }
}
