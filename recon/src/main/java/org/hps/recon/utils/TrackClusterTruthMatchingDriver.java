package org.hps.recon.utils;

import hep.aida.IAnalysisFactory;
import hep.aida.IHistogram1D;
import hep.aida.IHistogram2D;
import hep.aida.IHistogramFactory;
import hep.aida.ITree;
import hep.aida.ref.rootwriter.RootFileStore;

import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Vector;

import java.util.HashMap;
import java.util.List; 
import java.util.ArrayList; 
import java.util.Map;
import java.util.Set;
import java.util.HashSet;
import java.io.IOException;
import java.util.logging.Logger;

import org.hps.recon.tracking.TrackUtils;
import org.hps.recon.tracking.TrackData;
import org.hps.record.StandardCuts;
import org.hps.recon.ecal.cluster.ClusterUtilities;
//import org.hps.analysis.MC.TrackTruthMatching;

import org.hps.util.Pair;
import org.hps.util.RK4integrator;

import org.lcsim.event.EventHeader;
import org.lcsim.event.LCRelation;
import org.lcsim.event.MCParticle;
import org.lcsim.event.SimTrackerHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.Track;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.event.RelationalTable;
import org.lcsim.event.base.BaseRelationalTable;
import org.lcsim.util.Driver;
import org.lcsim.geometry.Detector;
import org.lcsim.event.TrackState;
import org.lcsim.event.SimCalorimeterHit;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.geometry.subdetector.HPSEcal3;

/** 
 * Driver used to truth match Tracks and Clusters, and then check the
 * performance of Track-to-Cluster Matcher algorithm.
 **/

public class TrackClusterTruthMatchingDriver extends Driver {

    //logger
    static private final Logger LOGGER = Logger.getLogger(TrackClusterMatcherNSigma.class.getPackage().getName());

    //fieldmap
    private org.lcsim.geometry.FieldMap fM;
    //Plot utilities
    private ITree tree;
    private IHistogramFactory histogramFactory;
    private Map<String, IHistogram1D> plots1D;
    private Map<String, IHistogram2D> plots2D;

    //enable truth plots
    boolean enablePlots = false;

    //Relational Tables
    RelationalTable hitToRotated = null;
    RelationalTable hitToStrips = null;
    RelationalTable TrktoData = new BaseRelationalTable(RelationalTable.Mode.ONE_TO_ONE, RelationalTable.Weighting.UNWEIGHTED);

    //inputs required for track cluster matching
    HPSEcal3 ecal;
    double beamEnergy = 2.3;
    protected StandardCuts cuts = new StandardCuts();

    //matcher
    TrackClusterMatcher matcher;
    private String clusterParamFileName = null;
    protected boolean enableTrackClusterMatchPlots = false;

    //misc
    private int flipSign = 1;
    protected double bField;

    //time difference for track and cluster
    double trackClusterTimeOffset;

    //scoring plane hits
    String ecalScoringPlaneHitsCollectionName = "TrackerHitsECal";
    String trackToScoringPlaneHitRelationsName = "TrackToEcalScoringPlaneHitRelations";
    //reco tracks collection name
    String trackCollectionName = "KalmanFullTracks";
    //ecal clusters
    String ecalClustersCollectionName = "EcalClusters";
    String ecalReadoutHitsCollectionName = "EcalReadoutHits";
    String ecalTruthRelationsName = "EcalTruthRelations";

    //matcher
    String trackClusterMatcherAlgo = "TrackClusterMatcherMinDistance";

    private Set<SimTrackerHit> simhitsontrack = new HashSet<SimTrackerHit>();

    public void setEcalClustersCollectionName(String name){
        this.ecalClustersCollectionName = name;
    }

    public void setTrackClusterTimeOffset(double input) {
        trackClusterTimeOffset = input;
        System.out.println("[truth] setting offset");
        cuts.setTrackClusterTimeOffset(input);
    }
    
    public void setTrackCollectionName(String trackCollectionName) {
        this.trackCollectionName = trackCollectionName;
    }
    
    public void setTrackClusterMatcherAlgo(String algoName) {
        this.trackClusterMatcherAlgo = algoName;
    }

    public void setEnableTrackClusterMatchPlots(boolean input){
        enableTrackClusterMatchPlots = input;
    }
    
    public void setBeamEnergy(double beamEnergy){
        this.beamEnergy = beamEnergy;
    }
    
    public void setClusterParamFileName(String input){
        this.clusterParamFileName = input;
    }

    public void setEnablePlots(boolean enablePlots){
        this.enablePlots = enablePlots;
    }

    public void setUseCorrectedClusterPositionsForMatching(boolean val) {
        useCorrectedClusterPositionsForMatching = val;
    }

    public void setUseTrackPositionForClusterCorrection(boolean val) {
        useTrackPositionForClusterCorrection = val;
    }

    public void setApplyClusterCorrections(boolean val) {
        applyClusterCorrections = val;
    }

    public void setMaxMatchDt(double input){
        cuts.setMaxMatchDt(input);
    }

    boolean useCorrectedClusterPositionsForMatching = false;
    boolean useTrackPositionForClusterCorrection = true;
    boolean applyClusterCorrections = true;

    /**
     * Updates the magnetic field parameters to match the appropriate values for
     * the current detector settings.
     */
    protected void detectorChanged(Detector detector) {

        if (clusterParamFileName == null) {
            if (beamEnergy > 2) {
                setClusterParamFileName("ClusterParameterization2016.dat");
            } else {
                setClusterParamFileName("ClusterParameterization2015.dat");
            }
        }

        //By default, use the original track-cluster matching class
        matcher = TrackClusterMatcherFactory.create(trackClusterMatcherAlgo);
        //matcher.initializeParameterization(clusterParamFileName);
        matcher.setBFieldMap(detector.getFieldMap());
        matcher.setTrackCollectionName(trackCollectionName);
        matcher.enablePlots(enableTrackClusterMatchPlots);

        // Set the magnetic field parameters to the appropriate values.
        Hep3Vector ip = new BasicHep3Vector(0., 0., 500.0);
        bField = detector.getFieldMap().getField(ip).y();
        if (bField < 0) {
            flipSign = -1;
        }

        ecal = (HPSEcal3) detector.getSubdetector("Ecal");

        if (cuts == null) {
            cuts = new StandardCuts(beamEnergy);
        } else {
            cuts.changeBeamEnergy(beamEnergy);
        }
    }

    public void saveHistograms() {
        String rootFile = String.format("%s_truthTrackClusterMatching.root",this.trackCollectionName);
        RootFileStore store = new RootFileStore(rootFile);
        try {
            store.open();
            store.add(tree);
            store.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void bookHistograms() {

        String trackCollectionName = this.trackCollectionName;
        plots1D = new HashMap<String, IHistogram1D>();
        plots2D = new HashMap<String, IHistogram2D>();
        tree = IAnalysisFactory.create().createTreeFactory().create();
        histogramFactory = IAnalysisFactory.create().createHistogramFactory(tree);

//Plots for checking new Track MCP matching tools
        plots1D.put(String.format("track_max_mcp_hit_multiplicity",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("track_max_mcp_hit_multiplicity",this.trackCollectionName), 100, 0, 100));
        plots1D.put(String.format("track_most-next_mcp_hit_multiplicity",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("track_most-next_mcp_hit_multiplicity",this.trackCollectionName), 100, 0, 100));
        plots1D.put(String.format("track_number_of_mcps",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("track_number_of_mcps",this.trackCollectionName), 40, 0, 40));
        plots2D.put(String.format("track_max_mcp_vs_next_best_mcp_hits",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("track_max_mcp_vs_next_best_mcp_hits",this.trackCollectionName), 40, 0, 40, 40, 0, 40));
        plots2D.put(String.format("track_max_mcp_vs_next_best_mcp_hits_same_pdgid",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("track_max_mcp_vs_next_best_mcp_hits_same_pdgid",this.trackCollectionName), 40, 0, 40, 40, 0, 40));
        plots2D.put(String.format("track_max_mcp_vs_next_best_mcp_hits_different_pdgid",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("track_max_mcp_vs_next_best_mcp_hits_different_pdgid",this.trackCollectionName), 40, 0, 40, 40, 0, 40));

        plots1D.put(String.format("number_strip_hits_per_sensor_layer",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("number_strip_hits_per_sensor_layer",this.trackCollectionName), 10, 0, 10));
        plots1D.put(String.format("number_of_mcps_on_striphits",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("number_of_mcps_on_striphits",this.trackCollectionName), 20, 0, 20));

        plots1D.put(String.format("best_mcp_nStripHits_over_total_nStripHits_on_track",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("best_mcp_nStripHits_over_total_nStripHits_on_track",this.trackCollectionName), 100, 0, 2));

        plots2D.put(String.format("best_mcp_nSensorsHit_v_2nd_best_mcp_nSensorsHit",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("best_mcp_nSensorsHit_v_2nd_best_mcp_nSensorsHit",this.trackCollectionName), 14, 0, 14, 14, 0, 14));

        plots2D.put(String.format("number_mcps_on_si_cluster",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("number_mcps_on_si_cluster",this.trackCollectionName), 13, 0, 13, 10, 0, 10));

        plots2D.put(String.format("new_trackMCP_match_trackP_v_mcpP",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("new_trackMCP_match_trackP_v_mcpP",this.trackCollectionName), 1000, -5, 5, 1000, -5, 5));
        plots2D.put(String.format("existing_trackMCP_match_trackP_v_mcpP",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("existing_trackMCP_match_trackP_v_mcpP",this.trackCollectionName), 1000, -5, 5, 1000, -5, 5));

//Plots for creating track+cluster residual parameterization files
        plots2D.put(String.format("%s_ele_TOP_track_cluster_truth_pairs_residual_param_dx",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_TOP_track_cluster_truth_pairs_residual_param_dx",this.trackCollectionName),50, 0, 5, 160,-40,40));
        plots2D.put(String.format("%s_ele_TOP_track_cluster_truth_pairs_residual_param_dy",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_TOP_track_cluster_truth_pairs_residual_param_dy",this.trackCollectionName),50, 0, 5, 160,-40,40));
        plots2D.put(String.format("%s_ele_TOP_track_cluster_truth_pairs_residual_param_dz",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_TOP_track_cluster_truth_pairs_residual_param_dz",this.trackCollectionName),50, 0, 5, 160,-40,40));
        plots2D.put(String.format("%s_ele_BOTTOM_track_cluster_truth_pairs_residual_param_dx",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_BOTTOM_track_cluster_truth_pairs_residual_param_dx",this.trackCollectionName),50, 0, 5, 160,-40,40));
        plots2D.put(String.format("%s_ele_BOTTOM_track_cluster_truth_pairs_residual_param_dy",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_BOTTOM_track_cluster_truth_pairs_residual_param_dy",this.trackCollectionName),50, 0, 5, 160,-40,40));
        plots2D.put(String.format("%s_ele_BOTTOM_track_cluster_truth_pairs_residual_param_dz",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_BOTTOM_track_cluster_truth_pairs_residual_param_dz",this.trackCollectionName),50, 0, 5, 160,-40,40));

        plots2D.put(String.format("%s_pos_TOP_track_cluster_truth_pairs_residual_param_dx",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_TOP_track_cluster_truth_pairs_residual_param_dx",this.trackCollectionName),50, 0, 5, 160,-40,40));
        plots2D.put(String.format("%s_pos_TOP_track_cluster_truth_pairs_residual_param_dy",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_TOP_track_cluster_truth_pairs_residual_param_dy",this.trackCollectionName),50, 0, 5, 160,-40,40));
        plots2D.put(String.format("%s_pos_TOP_track_cluster_truth_pairs_residual_param_dz",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_TOP_track_cluster_truth_pairs_residual_param_dz",this.trackCollectionName),50, 0, 5, 160,-40,40));
        plots2D.put(String.format("%s_pos_BOTTOM_track_cluster_truth_pairs_residual_param_dx",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_BOTTOM_track_cluster_truth_pairs_residual_param_dx",this.trackCollectionName),50, 0, 5, 160,-40,40));
        plots2D.put(String.format("%s_pos_BOTTOM_track_cluster_truth_pairs_residual_param_dy",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_BOTTOM_track_cluster_truth_pairs_residual_param_dy",this.trackCollectionName),50, 0, 5, 160,-40,40));
        plots2D.put(String.format("%s_pos_BOTTOM_track_cluster_truth_pairs_residual_param_dz",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_BOTTOM_track_cluster_truth_pairs_residual_param_dz",this.trackCollectionName),50, 0, 5, 160,-40,40));
    
//goodMatch plots

        plots2D.put(String.format("%s_ele_TOP_track_cluster_good_matches_dx",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_TOP_track_cluster_good_matches_dx",this.trackCollectionName),50, 0, 5, 160,-40,40));
        plots2D.put(String.format("%s_ele_TOP_track_cluster_good_matches_dy",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_TOP_track_cluster_good_matches_dy",this.trackCollectionName),50, 0, 5, 160,-40,40));
        plots2D.put(String.format("%s_ele_BOTTOM_track_cluster_good_matches_dx",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_BOTTOM_track_cluster_good_matches_dx",this.trackCollectionName),50, 0, 5, 160,-40,40));
        plots2D.put(String.format("%s_ele_BOTTOM_track_cluster_good_matches_dy",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_BOTTOM_track_cluster_good_matches_dy",this.trackCollectionName),50, 0, 5, 160,-40,40));
        plots2D.put(String.format("%s_ele_good_matches_trackP_v_clusterE",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_good_matches_trackP_v_clusterE",this.trackCollectionName),500, 0, 5, 500, 0, 5));

        plots2D.put(String.format("%s_pos_TOP_track_cluster_good_matches_dx",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_TOP_track_cluster_good_matches_dx",this.trackCollectionName),50, 0, 5, 160,-40,40));
        plots2D.put(String.format("%s_pos_TOP_track_cluster_good_matches_dy",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_TOP_track_cluster_good_matches_dy",this.trackCollectionName),50, 0, 5, 160,-40,40));
        plots2D.put(String.format("%s_pos_BOTTOM_track_cluster_good_matches_dx",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_BOTTOM_track_cluster_good_matches_dx",this.trackCollectionName),50, 0, 5, 160,-40,40));
        plots2D.put(String.format("%s_pos_BOTTOM_track_cluster_good_matches_dy",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_BOTTOM_track_cluster_good_matches_dy",this.trackCollectionName),50, 0, 5, 160,-40,40));
        plots2D.put(String.format("%s_pos_good_matches_trackP_v_clusterE",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_good_matches_trackP_v_clusterE",this.trackCollectionName),500, 0, 5, 500, 0, 5));

//Total Counts of Interesting Events
    
        plots1D.put(String.format("nMCPsEvaluated",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nMCPsEvaluated",this.trackCollectionName), 101, -1, 100));
        plots1D.put(String.format("nGoodMatches",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nGoodMatches",this.trackCollectionName), 101, -1, 100));
        plots1D.put(String.format("nBadMatches",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nBadMatches",this.trackCollectionName), 101, -1, 100));
        plots1D.put(String.format("nMissedMatches",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nMissedMatches",this.trackCollectionName), 101, -1, 100));
        plots1D.put(String.format("nUnmatchedAmbiguousClusterInefficiency",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nUnmatchedAmbiguousClusterInefficiency",this.trackCollectionName), 101, -1, 100));
        plots1D.put(String.format("nMatchedAmbiguousTrackInefficiency",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nMatchedAmbiguousTrackInefficiency",this.trackCollectionName), 101, -1, 100));
        plots1D.put(String.format("nMatchedAmbiguousClusterInefficiency",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nMatchedAmbiguousClusterInefficiency",this.trackCollectionName), 101, -1, 100));
        plots1D.put(String.format("nUnmatchedAmbiguousTrackInefficiency",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nUnmatchedAmbiguousTrackInefficiency",this.trackCollectionName), 101, -1, 100));
        plots1D.put(String.format("nPhotons_matched_to_track",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nPhotons_matched_to_track",this.trackCollectionName), 101, -1, 100));
        plots1D.put(String.format("nTruth_photon_clusters",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nTruth_photon_clusters",this.trackCollectionName), 101, -1, 100));
        plots1D.put(String.format("nRogueMatches",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nRogueMatches",this.trackCollectionName), 101, -1, 100));
        plots1D.put(String.format("nRogueTracks",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nRogueTracks",this.trackCollectionName), 101, -1, 100));
        plots1D.put(String.format("nIddPhotons",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nIddPhotons",this.trackCollectionName), 101, -1, 100));

        plots1D.put(String.format("nTracklessMCPClusterMatchedToRogueTrack",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nTracklessMCPClusterMatchedToRogueTrack",this.trackCollectionName), 101, -1, 100));
        plots1D.put(String.format("nTracklessMCPCorrectlyUnmatched",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nTracklessMCPCorrectlyUnmatched",this.trackCollectionName), 101, -1, 100));
        plots1D.put(String.format("nClusterlessMCPMatchedToRogueCluster",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nClusterlessMCPMatchedToRogueCluster",this.trackCollectionName), 101, -1, 100));
        plots1D.put(String.format("nClusterlessMCPCorrectlyUnmatched",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nClusterlessMCPCorrectlyUnmatched",this.trackCollectionName), 101, -1, 100));
        plots1D.put(String.format("nMCParticlesEvaluated",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nMCParticlesEvaluated",this.trackCollectionName), 101, -1, 100));

//Event characterization plots

        plots1D.put(String.format("nEvents",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nEvents",this.trackCollectionName),10000 , 0, 10000));

        plots1D.put(String.format("nMCParticles_per_event",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nMCParticles_per_event",this.trackCollectionName),1000 , 0, 10000));

        plots1D.put(String.format("nClusters_per_event",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nClusters_per_event",this.trackCollectionName),100 , 0, 100));

        plots1D.put(String.format("nTracks_per_event",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nTracks_per_event",this.trackCollectionName),100 , 0, 100));

        plots2D.put(String.format("nTracks_v_nClusters_per_event",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("nTracks_v_nClusters_per_event",this.trackCollectionName),100, 0, 100, 100, 0, 100));


        plots2D.put(String.format("nMCP_w_truth_tracks_v_nMCP_w_truth_clusters_per_event",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("nMCP_w_truth_tracks_v_nMCP_w_truth_clusters_per_event",this.trackCollectionName),100, 0, 100, 100, 0, 100));
        plots1D.put(String.format("mcp_w_truth_cluster_but_no_truth_tracks_momentum",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcp_w_truth_cluster_but_no_truth_tracks_momentum",this.trackCollectionName),500, 0, 5));

        plots1D.put(String.format("mcp_w_truth_cluster_but_no_truth_tracks_nSimTrackerHits",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcp_w_truth_cluster_but_no_truth_tracks_nSimTrackerHits",this.trackCollectionName),30, 0, 30));

        plots1D.put(String.format("mcp_w_truth_tracks_but_no_truth_cluster_momentum",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcp_w_truth_tracks_but_no_truth_cluster_momentum",this.trackCollectionName),500, 0, 5));

        plots1D.put(String.format("mcp_w_truth_cluster_but_no_truth_tracks_nSimCalHits",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcp_w_truth_cluster_but_no_truth_tracks_nSimCalHits",this.trackCollectionName), 30, 0, 30));

        plots1D.put(String.format("mcp_w_truth_tracks_AND_truth_cluster_momentum",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcp_w_truth_tracks_AND_truth_cluster_momentum",this.trackCollectionName),500, 0, 5));

        plots2D.put(String.format("mcp_w_truth_track_AND_truth_cluster_track_v_cluster_momentum",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("mcp_w_truth_track_AND_truth_cluster_track_v_cluster_momentum",this.trackCollectionName),500, 0, 5, 500, 0, 5));

        plots1D.put(String.format("mcps_all_checked_against_matcher_momentum",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcps_all_checked_against_matcher_momentum",this.trackCollectionName),500, 0, 5));

        plots2D.put(String.format("mcps_loop_nSimTrackerHits_v_nSimCalHits",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("mcps_loop_nSimTrackerHits_v_nSimCalHits",this.trackCollectionName),30, 0, 30, 30, 0, 30));

        plots1D.put(String.format("nMCP_w_truth_tracks_BUT_NO_truth_cluster_per_event",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nMCP_w_truth_tracks_BUT_NO_truth_cluster_per_event",this.trackCollectionName),100 , 0, 100));

        plots1D.put(String.format("nMCP_w_truth_cluster_BUT_NO_truth_tracks_per_event",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nMCP_w_truth_cluster_BUT_NO_truth_tracks_per_event",this.trackCollectionName),100 , 0, 100));

        plots1D.put(String.format("nMCP_w_truth_tracks_AND_truth_cluster_per_event",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nMCP_w_truth_tracks_AND_truth_cluster_per_event",this.trackCollectionName),100 , 0, 100));


//MCP LOOP OVER MATCHE RESULTS

        plots2D.put(String.format("mcp_loop_truth_cluster_matched_to_wrong_track_p_v_p",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("mcp_loop_truth_cluster_matched_to_wrong_track_p_v_p",this.trackCollectionName),500, 0, 5, 500, 0, 5));

        plots1D.put(String.format("mcp_loop_trackless_truth_cluster_matched_to_rogue_track_nSimTrackerHits",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcp_loop_trackless_truth_cluster_matched_to_rogue_track_nSimTrackerHits",this.trackCollectionName),30 , 0, 30));

        plots1D.put(String.format("mcp_loop_no_truth_tracks_truth_cluster_not_matched_to_track_nSimTrackerHits",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcp_loop_no_truth_tracks_truth_cluster_not_matched_to_track_nSimTrackerHits",this.trackCollectionName),30 , 0, 30));

        plots2D.put(String.format("mcp_loop_truth_track_matched_to_wrong_cluster_E_v_E",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("mcp_loop_truth_track_matched_to_wrong_cluster_E_v_E",this.trackCollectionName),500, 0, 5, 500, 0, 5));



        plots1D.put(String.format("mcp_loop_nEle_per_event",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcp_loop_nEle_per_event",this.trackCollectionName),100 , 0, 100));
        
        plots1D.put(String.format("mcp_loop_nPos_per_event",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcp_loop_nPos_per_event",this.trackCollectionName),100 , 0, 100));
        
        plots1D.put(String.format("mcp_loop_nPhotons_per_event",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcp_loop_nPhotons_per_event",this.trackCollectionName),100 , 0, 100));
        
        plots1D.put(String.format("mcp_loop_nCharged_particles_per_event",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcp_loop_nCharged_particles_per_event",this.trackCollectionName),100 , 0, 100));
        
        plots2D.put(String.format("mcp_loop_nCharged_particles_v_nTracks_per_event",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("mcp_loop_nCharged_particles_v_nTracks_per_event",this.trackCollectionName),100, 0, 100, 100, 0, 100));
        
        plots2D.put(String.format("mcp_loop_nMCPs_v_nClusters_per_event",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("mcp_loop_nMCPs_v_nClusters_per_event",this.trackCollectionName),100, 0, 100, 100, 0, 100));

//MATCHING PLOT

        plots1D.put(String.format("%s_mcp_truth_track_cluster_pair_momentum",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_mcp_truth_track_cluster_pair_momentum",this.trackCollectionName), 500, 0, 5));
        plots1D.put(String.format("%s_mcp_truth_track_cluster_pair_energy",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_mcp_truth_track_cluster_pair_energy",this.trackCollectionName), 500, 0, 5));

        plots2D.put(String.format("mcp_truth_track_cluster_pair_ntrackersimhits_v_momentum",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("mcp_truth_track_cluster_pair_ntrackersimhits_v_momentum",this.trackCollectionName),20, 0, 20, 500, 0, 5));
        plots2D.put(String.format("mcp_truth_track_ntrackersimhits_v_mcp_momentum",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("mcp_truth_track_ntrackersimhits_v_mcp_momentum",this.trackCollectionName), 20, 0, 20, 500, 0, 5));

//MCP TRUTH MATCH TO TRACKS

        plots1D.put(String.format("mcp_nHits_in_Tracker",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcp_nHits_in_Tracker",this.trackCollectionName), 20, 0, 20));
        
        plots1D.put(String.format("mcp_ofInterest_nRawTrackerHits_per_track",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcp_ofInterest_nRawTrackerHits_per_track",this.trackCollectionName), 30, 0, 30));

        plots2D.put(String.format("mcp_ofInterest_momentum_v_nRawTrackerHits_per_track",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("mcp_ofInterest_momentum_v_nRawTrackerHits_per_track",this.trackCollectionName),500, 0, 5, 30, 0, 30));

        plots1D.put(String.format("mcp_mostHitsTrack_nHits",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcp_mostHitsTrack_nHits",this.trackCollectionName), 20, 0, 20));

        plots1D.put(String.format("mcp_nTracks",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcp_nTracks",this.trackCollectionName), 20, 0, 20));

        plots2D.put(String.format("mcp_momentum_v_nTracks",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("mcp_momentum_v_nTracks",this.trackCollectionName),500, 0, 5, 20, 0, 20));

        plots1D.put(String.format("nMCP_primary_FEEs_per_event",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nMCP_primary_FEEs_per_event",this.trackCollectionName), 500, 0, 500));

        plots1D.put(String.format("nMCP_radiative_FEEs_per_event",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nMCP_radiative_FEEs_per_event",this.trackCollectionName), 500, 0, 500));

        plots1D.put(String.format("nMCP_622_primary_daughters_per_event",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nMCP_622_primary_daughters_per_event",this.trackCollectionName), 500, 0, 500));

        plots1D.put(String.format("nMCP_623_primary_daughters_per_event",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nMCP_623_primary_daughters_per_event",this.trackCollectionName), 500, 0, 500));

        plots2D.put(String.format("nMCP_622or623_primary_daughters_v_primary_FEEs_per_event",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("nMCP_622or623_primary_daughters_v_primary_FEEs_per_event",this.trackCollectionName),500, 0, 500, 500, 0, 500));

        plots2D.put(String.format("nMCP_622and623_primary_daughters_and_primary_FEEs_v_nTracks_per_event",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("nMCP_622and623_primary_daughters_and_primary_FEEs_v_nTracks_per_event",this.trackCollectionName),500, 0, 500, 500, 0, 500));

        plots1D.put(String.format("nMCPs_on_track",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("nMCPs_on_track",this.trackCollectionName), 20, 0, 20));

        plots1D.put(String.format("tracks_unmatched_to_mcp_ofInterest_momentum",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("tracks_unmatched_to_mcp_ofInterest_momentum",this.trackCollectionName), 500, 0, 5));

        plots2D.put(String.format("mcp_nSimTrackerHits_v_sum_nHitsOnTrack",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("mcp_nSimTrackerHits_v_sum_nHitsOnTrack",this.trackCollectionName),20, 0, 20, 20, 0, 20));

        plots1D.put(String.format("mcp_FEE_nTracks_disambiguated",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcp_FEE_nTracks_disambiguated",this.trackCollectionName), 20, 0, 20));

        plots1D.put(String.format("mcp_nonFEE_nTracks_disambiguated",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcp_nonFEE_nTracks_disambiguated",this.trackCollectionName), 20, 0, 20));

        plots1D.put(String.format("mcp_nTracks_disambiguated",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcp_nTracks_disambiguated",this.trackCollectionName), 20, 0, 20));

        plots2D.put(String.format("mcp_momentum_v_nTracks_disambiguated",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("mcp_momentum_v_nTracks_disambiguated",this.trackCollectionName),500, 0, 5, 20, 0, 20));

        plots1D.put(String.format("mcp_nhits_on_track_disambiguated",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcp_nhits_on_track_disambiguated",this.trackCollectionName), 20, 0, 20));

        plots1D.put(String.format("mcp_most_hits_on_track_disambiguated",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcp_most_hits_on_track_disambiguated",this.trackCollectionName), 20, 0, 20));

        plots1D.put(String.format("mcp_FEEs_most_hits_on_track_disambiguated",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcp_FEEs_most_hits_on_track_disambiguated",this.trackCollectionName), 20, 0, 20));

        plots1D.put(String.format("mcp_NOT_FEEs_most_hits_on_track_disambiguated",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcp_NOT_FEEs_most_hits_on_track_disambiguated",this.trackCollectionName), 20, 0, 20));

        plots2D.put(String.format("mcp_momentum_v_best_track_momentum_disam",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("mcp_momentum_v_best_track_momentum_disam",this.trackCollectionName),500, 0, 5, 500, 0, 5));

        plots1D.put(String.format("mcp_momentum_best_track_momentum_ratio_disam",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcp_momentum_best_track_momentum_ratio_disam",this.trackCollectionName), 200, 0, 2));

        plots2D.put(String.format("mcp_nSimTrackerHits_most_nHitsOnTrack_disam",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("mcp_nSimTrackerHits_most_nHitsOnTrack_disam",this.trackCollectionName),20, 0, 20, 20, 0, 20));

        plots1D.put(String.format("mcp_FEE_wAtLeast_one_track_per_event_disamb",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcp_FEE_wAtLeast_one_track_per_event_disamb",this.trackCollectionName), 500, 0, 500));

        plots1D.put(String.format("mcp_NOT_FEE_wAtLeast_one_track_per_event_disamb",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("mcp_NOT_FEE_wAtLeast_one_track_per_event_disamb",this.trackCollectionName), 500, 0, 500));




//FINAL CLUSTER TRUTH MATCHING PLOTS

        //Cluster Positions in XY plane
        plots2D.put(String.format("ecal_cluster_positions_xy_plane",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("ecal_cluster_positions_xy_plane",this.trackCollectionName),1000, -500, 500,1000, -500, 500));

        plots2D.put(String.format("cluster_truthpositions_xy_plane",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("cluster_truthpositions_xy_plane",this.trackCollectionName),50, -280, 370, 18, -110, 124));
        
        plots1D.put(String.format("cluster_truth_stage_0_energy",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_truth_stage_0_energy",this.trackCollectionName), 1000, 0, 10));

        plots1D.put(String.format("cluster_truth_stage_1_energy",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_truth_stage_1_energy",this.trackCollectionName), 1000, 0, 10));

        plots2D.put(String.format( "cluster_truth_stage_1_mcpEndpointz_v_ds"), histogramFactory.createHistogram2D(String.format("cluster_truth_stage_1_mcpEndpointz_v_ds"),200, -100, 1900,1000, 0, 2000));

        plots1D.put(String.format("cluster_truth_stage_1_energy_ratio",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_truth_stage_1_energy_ratio",this.trackCollectionName), 1000, 0, 10));
        
        plots2D.put(String.format("cluster_truth_stage_1_cluster_v_mcp_energy",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("cluster_truth_stage_1_cluster_v_mcp_energy",this.trackCollectionName),1000, 0, 5,1000, 0, 5));

        plots1D.put(String.format("cluster_truth_stage_2_energy",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_truth_stage_2_energy",this.trackCollectionName), 1000, 0, 10));
        
        plots2D.put(String.format( "cluster_truth_stage_2_mcpEndpointz_v_ds"), histogramFactory.createHistogram2D(String.format("cluster_truth_stage_2_mcpEndpointz_v_ds"),200, -100, 1900,1000, 0, 2000));
        
        plots1D.put(String.format("cluster_truth_stage_2_energy_ratio",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_truth_stage_2_energy_ratio",this.trackCollectionName), 1000, 0, 10));
        
        plots2D.put(String.format("cluster_truth_stage_2_cluster_v_mcp_energy",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("cluster_truth_stage_2_cluster_v_mcp_energy",this.trackCollectionName),1000, 0, 5,1000, 0, 5));

        plots2D.put(String.format("cluster_truth_stage_2_mcpPy_v_clustery",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("cluster_truth_stage_2_mcpPy_v_clustery",this.trackCollectionName),1000, -5, 5, 1000, -500, 500));

        plots1D.put(String.format("cluster_truth_stage_3_energy",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_truth_stage_3_energy",this.trackCollectionName), 1000, 0, 10));
        
        plots2D.put(String.format( "cluster_truth_stage_3_mcpEndpointz_v_ds"), histogramFactory.createHistogram2D(String.format("cluster_truth_stage_3_mcpEndpointz_v_ds"),200, -100, 1900,1000, 0, 2000));
        
        plots1D.put(String.format("cluster_truth_stage_3_energy_ratio",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_truth_stage_3_energy_ratio",this.trackCollectionName), 1000, 0, 10));
        
        plots2D.put(String.format("cluster_truth_stage_3_cluster_v_mcp_energy",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("cluster_truth_stage_3_cluster_v_mcp_energy",this.trackCollectionName),1000, 0, 5,1000, 0, 5));

        plots1D.put(String.format("clusters_matched_to_n_mcps"), histogramFactory.createHistogram1D(String.format("clusters_matched_to_n_mcps"),  10, 0, 10));

        //Cluster Duplicates
        plots1D.put(String.format("cluster_truth_stage_3_duplicate_mcp_match_dx",trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_truth_stage_3_duplicate_mcp_match_dx",trackCollectionName),  1000, -1000, 1000));
        
        plots1D.put(String.format("cluster_truth_stage_3_duplicate_mcp_match_dy",trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_truth_stage_3_duplicate_mcp_match_dy",trackCollectionName),  1000, -1000, 1000));

        plots1D.put(String.format("cluster_truth_stage_3_cut_remaining_duplicate_mcp_match_dx",trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_truth_stage_3_cut_remaining_duplicate_mcp_match_dx",trackCollectionName),  1000, -1000, 1000));
        
        plots1D.put(String.format("cluster_truth_stage_3_cut_remaining_duplicate_mcp_match_dy",trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_truth_stage_3_cut_remaining_duplicate_mcp_match_dy",trackCollectionName),  1000, -1000, 1000));

        //Final Cut Clusters
        plots1D.put(String.format("cluster_truth_stage_final_energy",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_truth_stage_final_energy",this.trackCollectionName), 1000, 0, 10));
        
        plots2D.put(String.format( "cluster_truth_stage_final_mcpEndpointz_v_ds"), histogramFactory.createHistogram2D(String.format("cluster_truth_stage_final_mcpEndpointz_v_ds"),200, -100, 1900,1000, 0, 2000));
        
        plots1D.put(String.format("cluster_truth_stage_final_energy_ratio",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_truth_stage_final_energy_ratio",this.trackCollectionName), 1000, 0, 10));
        
        plots2D.put(String.format("cluster_truth_stage_final_cluster_v_mcp_energy",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("cluster_truth_stage_final_cluster_v_mcp_energy",this.trackCollectionName),1000, 0, 5,1000, 0, 5));

        plots1D.put(String.format("cluster_truth_energy",trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_truth_energy",trackCollectionName),  100, 0, 5));



        //Plot the XY acceptance of the Ecal by adding half the crystal width
        //to the cluster xy positions
        plots2D.put(String.format("ecal_crystal_acceptance_xy"), histogramFactory.createHistogram2D(String.format("ecal_crystal_acceptance_xy"),800, -400, 400, 240, -120, 120));


//SCORING PLANE PLOTS
        //Plots showing residuals between track at Ecal and truth-matched
        //scoring plane hit that has been extrapolated to the Ecal
        plots1D.put(String.format("%s_ele_track_scoringplane_hit_dx",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_track_scoringplane_hit_dx",this.trackCollectionName), 800, -400, 400));
        plots1D.put(String.format("%s_ele_track_scoringplane_hit_dy",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_track_scoringplane_hit_dy",this.trackCollectionName), 800, -400, 400));
        plots1D.put(String.format("%s_ele_track_scoringplane_hit_dz",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_track_scoringplane_hit_dz",this.trackCollectionName), 800, -400, 400));
        plots1D.put(String.format("%s_ele_track_scoringplane_hit_dr",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_track_scoringplane_hit_dr",this.trackCollectionName), 800, -400, 400));
        plots1D.put(String.format("%s_ele_track_scoringplane_hit_dt",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_track_scoringplane_hit_dt",this.trackCollectionName), 800, -400, 400));

        plots1D.put(String.format("%s_pos_track_scoringplane_hit_dx",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_track_scoringplane_hit_dx",this.trackCollectionName), 800, -400, 400));
        plots1D.put(String.format("%s_pos_track_scoringplane_hit_dy",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_track_scoringplane_hit_dy",this.trackCollectionName), 800, -400, 400));
        plots1D.put(String.format("%s_pos_track_scoringplane_hit_dz",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_track_scoringplane_hit_dz",this.trackCollectionName), 800, -400, 400));
        plots1D.put(String.format("%s_pos_track_scoringplane_hit_dr",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_track_scoringplane_hit_dr",this.trackCollectionName), 800, -400, 400));
        plots1D.put(String.format("%s_pos_track_scoringplane_hit_dt",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_track_scoringplane_hit_dt",this.trackCollectionName), 800, -400, 400));

            //scoringplane momentum components
        plots1D.put(String.format("%s_ele_track_scoringplane_hit_px",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_track_scoringplane_hit_px",this.trackCollectionName), 600, -1, 5));
        plots1D.put(String.format("%s_ele_track_scoringplane_hit_py",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_track_scoringplane_hit_py",this.trackCollectionName), 600, -1, 5));
        plots1D.put(String.format("%s_ele_track_scoringplane_hit_pz",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_track_scoringplane_hit_pz",this.trackCollectionName), 600, -1, 5));

        plots1D.put(String.format("%s_pos_track_scoringplane_hit_px",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_track_scoringplane_hit_px",this.trackCollectionName), 600, -1, 5));
        plots1D.put(String.format("%s_pos_track_scoringplane_hit_py",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_track_scoringplane_hit_py",this.trackCollectionName), 600, -1, 5));
        plots1D.put(String.format("%s_pos_track_scoringplane_hit_pz",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_track_scoringplane_hit_pz",this.trackCollectionName), 600, -1, 5));
            //track v scoringplane momentum
        plots2D.put(String.format("%s_ele_scoringplaneHit_v_truth_track_p",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_scoringplaneHit_v_truth_track_p",this.trackCollectionName),600, -1, 5, 600, -1, 5));
        plots2D.put(String.format("%s_pos_scoringplaneHit_v_truth_track_p",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_scoringplaneHit_v_truth_track_p",this.trackCollectionName),600, -1, 5, 600, -1, 5));


        //Track extrapolation to Ecal: Momentum vs truth-extrap position
        //residuals
        plots2D.put(String.format("%s_ele_RK4_scoringplanehit_to_ecal_ZvP",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_RK4_scoringplanehit_to_ecal_ZvP",this.trackCollectionName), 1500, 0, 1500,300,0,3));

        plots2D.put(String.format("%s_pos_RK4_scoringplanehit_to_ecal_ZvP",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_RK4_scoringplanehit_to_ecal_ZvP",this.trackCollectionName), 1500, 0, 1500,300,0,3));

//TRACK RECONSTRUCTION EFFICIENCY, ETC.

        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_dx",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_dx",trackCollectionName),  800, -200, 200));
        
        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_dy",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_dy",trackCollectionName),  800, -200, 200));
        
        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_dx",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_dx",trackCollectionName),  800, -200, 200));
        
        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_dy",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_dy",trackCollectionName),  800, -200, 200));
        
        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_dz",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_dz",trackCollectionName),  800, -200, 200));
        
        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_dz",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_dz",trackCollectionName),  800, -200, 200));
        
        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_EdivP",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_EdivP",trackCollectionName),  1000, 0, 10));
        
        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_EdivP",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_EdivP",trackCollectionName),  1000, 0, 10));

        plots2D.put(String.format("%s_truth_track_cluster_pair_xy_for_EdivP_gt_0.35_lt_0.6",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_truth_track_cluster_pair_xy_for_EdivP_gt_0.35_lt_0.6",this.trackCollectionName),800, -400, 400, 180, -90, 90));
        plots2D.put(String.format("%s_truth_track_cluster_pair_xy_for_EdivP_lt_0.35",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_truth_track_cluster_pair_xy_for_EdivP_lt_0.35",this.trackCollectionName),800, -400, 400, 180, -90, 90));
        plots2D.put(String.format("%s_truth_track_cluster_pair_xy_for_EdivP_gt_0.6",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_truth_track_cluster_pair_xy_for_EdivP_gt_0.6",this.trackCollectionName),800, -400, 400, 180, -90, 90));

        //Hit multiplicity for truth matching
        plots1D.put(String.format("%s_ele_track_maxMCPmultiplicity",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_track_maxMCPmultiplicity",this.trackCollectionName),16, 0, 16));
        plots1D.put(String.format("%s_pos_track_maxMCPmultiplicity",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_track_maxMCPmultiplicity",this.trackCollectionName),16, 0, 16));
        plots2D.put(String.format("%s_ele_track_p_v_maxMCPmultiplicity",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_track_p_v_maxMCPmultiplicity",this.trackCollectionName), 1000, 0, 5, 16, 0, 16));
        plots2D.put(String.format("%s_pos_track_p_v_maxMCPmultiplicity",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_track_p_v_maxMCPmultiplicity",this.trackCollectionName), 1000, 0, 5, 16, 0, 16));
        plots2D.put(String.format("%s_pos_track_p_v_mcp_p",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_track_p_v_mcp_p",this.trackCollectionName), 1000, 0, 5, 1000, 0, 5));
        plots2D.put(String.format("%s_ele_track_p_v_mcp_p",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_track_p_v_mcp_p",this.trackCollectionName), 1000, 0, 5, 1000, 0, 5));


        //track positions at Ecal in XY plane
        plots2D.put(String.format("%s_ele_track_xypos_at_ecal",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_track_xypos_at_ecal",this.trackCollectionName),1000, -500, 500,1000, -500, 500));
        plots2D.put(String.format("%s_pos_track_xypos_at_ecal",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_track_xypos_at_ecal",this.trackCollectionName),1000, -500, 500,1000, -500, 500));

        plots2D.put(String.format("%s_ele_truth_track_xypos_at_ecal",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_truth_track_xypos_at_ecal",this.trackCollectionName),1000, -500, 500,1000, -500, 500));
        plots2D.put(String.format("%s_pos_truth_track_xypos_at_ecal",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_truth_track_xypos_at_ecal",this.trackCollectionName),1000, -500, 500,1000, -500, 500));


        //Check track quality
        plots1D.put(String.format("%s_ele_track_chi2divndf",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_track_chi2divndf",this.trackCollectionName), 200, 0, 200));
        plots1D.put(String.format("%s_pos_track_chi2divndf",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_track_chi2divndf",this.trackCollectionName), 200, 0, 200));

    }

    public void startOfData() {

        bookHistograms();
    }

    public void endOfData() {

        if(enablePlots)
            matcher.saveHistograms();
        saveHistograms();
    }

    protected void process(EventHeader event) {

        //count objects event by event
        plots1D.get("nEvents").fill(1);
        
        //MCPs
        List<MCParticle> allmcps = event.get(MCParticle.class,"MCParticle");
        plots1D.get("nMCParticles_per_event").fill(allmcps.size());

        //Get EcalClusters from event
        List<Cluster> clusters = event.get(Cluster.class, ecalClustersCollectionName);
        plots1D.get("nClusters_per_event").fill(clusters.size());

        // Get collection of tracks from event
        List<Track> tracks = event.get(Track.class, trackCollectionName);
        plots1D.get("nTracks_per_event").fill(tracks.size());

        plots2D.get("nTracks_v_nClusters_per_event").fill(tracks.size(),clusters.size());

        //Apply time offset
        cuts.setTrackClusterTimeOffset(this.trackClusterTimeOffset);

        //Truth Match Clusters to MCParticles    
        Map<MCParticle, Cluster> mcpClustersMap = new HashMap<MCParticle, Cluster>();
        Map<Cluster, MCParticle> clustersMCPMap = new HashMap<Cluster, MCParticle>();
        getClusterMcpMap(clusters, event, false, mcpClustersMap, clustersMCPMap);

        //Plot cluster positions (move this to new class)
        List<Cluster> truthClusters = new ArrayList<Cluster>();
        for(Map.Entry<Cluster,MCParticle> entry : clustersMCPMap.entrySet()){
            Cluster truthcluster = entry.getKey();
            double clusterx = truthcluster.getPosition()[0];
            double clustery = truthcluster.getPosition()[1];
            double clusterz = truthcluster.getPosition()[2];
            double clusterEnergy = truthcluster.getEnergy();
            double energyRatio = clusterEnergy/clustersMCPMap.get(truthcluster).getEnergy();

            plots1D.get(String.format("cluster_truth_energy",trackCollectionName)).fill(clusterEnergy);
            plots2D.get(String.format("cluster_truthpositions_xy_plane")).fill(clusterx,clustery);
            truthClusters.add(entry.getKey());
        }

        /*
        drawEcalFace(truthClusters);

        //Get Kalman Track Data
        hitToRotated = TrackUtils.getHitToRotatedTable(event);
        hitToStrips = TrackUtils.getHitToStripsTable(event);
        List<TrackData> TrackData;
        List<LCRelation> trackRelations;
        if (this.trackCollectionName.contains("KalmanFullTracks")) {
            TrackData = event.get(TrackData.class, "KFTrackData");
            trackRelations = event.get(LCRelation.class, "KFTrackDataRelations");
            for (LCRelation relation : trackRelations) {
                if (relation != null && relation.getTo() != null){
                    TrktoData.add(relation.getFrom(), relation.getTo());
                }
            }
        }

        //Relational Table for getting Track truth information
        RelationalTable rawtomc = new BaseRelationalTable(RelationalTable.Mode.MANY_TO_MANY, RelationalTable.Weighting.UNWEIGHTED);
        if (event.hasCollection(LCRelation.class, "SVTTrueHitRelations")) {
            List<LCRelation> trueHitRelations = event.get(LCRelation.class, "SVTTrueHitRelations");
            for (LCRelation relation : trueHitRelations)
                if (relation != null && relation.getFrom() != null && relation.getTo() != null)
                    rawtomc.add(relation.getFrom(), relation.getTo());
        }

        //Initialize Truth Maps for MCParticles and Tracks.
        //Tracks are uniquely matched to one MCParticle.
        //MCParticles can be matched to multiple Tracks.
        Map<MCParticle, Map<Track, Integer>> mcpTrackMapDisamb = new HashMap<MCParticle, Map<Track, Integer>>();
        Map<Track, MCParticle> trackMCPMapDisamb = new HashMap<Track, MCParticle>();
        buildTrackTruthMaps(tracks, rawtomc, mcpTrackMapDisamb, trackMCPMapDisamb);

        //Init map that has full truth info for MCP
        Map<MCParticle, Pair<Cluster, List<Track>>> mcpFullTruthMap = new HashMap<MCParticle, Pair<Cluster, List<Track>>>();
        //Init alternative map that only keeps best Track, instead of
        //List<Track>
        Map<MCParticle, Pair<Track, Cluster>> mcpTrackClusterPairs = new HashMap<MCParticle, Pair<Track, Cluster>>();

        //Loop over truth Tracks
        for(Map.Entry<MCParticle, Map<Track, Integer>> entry : mcpTrackMapDisamb.entrySet()){
            MCParticle mcp = entry.getKey();
            List<Track> mcpTracks = new ArrayList<Track>();
            mcpTracks.addAll(entry.getValue().keySet());
            Track bestTrack = getMcpBestTrack(mcp,mcpTrackMapDisamb);
            int charge = -1* (int)Math.signum(bestTrack.getTrackStates().get(0).getOmega());
            double trackPmag = new BasicHep3Vector(bestTrack.getTrackStates().get(0).getMomentum()).magnitude();

            //If mcp has both truth Track(s) AND truth Cluster
            if(mcpClustersMap.containsKey(mcp)){
                Cluster mcpCluster = mcpClustersMap.get(mcp);
                Pair<Cluster, List<Track>> pair = new Pair<Cluster,List<Track>>(mcpCluster, mcpTracks);
                mcpFullTruthMap.put(mcp, pair);
                mcpTrackClusterPairs.put(mcp, new Pair<Track,Cluster>(bestTrack, mcpCluster));
                plots2D.get("mcp_matched_trackP_v_clusterE_x_mcp_charge").fill(charge*trackPmag, ((int)mcp.getCharge())*mcpCluster.getEnergy());
            }
            //If mcp has truth Track(s) but no truth Cluster is found
            else{
                plots2D.get("clusterless_mcpP_v_trackP").fill(((int)mcp.getCharge())*mcp.getMomentum().magnitude(), charge*trackPmag);   
                Pair<Cluster, List<Track>> pair = new Pair<Cluster,List<Track>>(null, mcpTracks);
                mcpFullTruthMap.put(mcp, pair);
            }
        }

        //Loop over truth Clusters that have no truth Tracks
        for(Map.Entry<MCParticle, Cluster> entry : mcpClustersMap.entrySet()){
            MCParticle mcp = entry.getKey();
            Cluster cluster = entry.getValue();
            //we've already added cases where both Track(s) and Cluster are
            //found. Skip.
            if(mcpTrackMapDisamb.containsKey(mcp))
                continue;

            Pair<Cluster, List<Track>> pair = new Pair<Cluster,List<Track>>(cluster, null);
            mcpFullTruthMap.put(mcp, pair);
            plots2D.get("trackless_mcpE_v_clusterE").fill(mcp.getEnergy(), cluster.getEnergy()); 
        }

        //Caluclate track+cluster pair position residuals as function of
        //momentum and create residual parameterization file to be used in
        //matcher algorithm class
        trackClusterResidualParameterization(mcpTrackClusterPairs);

        //Feed Track and Cluster collections to matcher algorithm to get
        //algorithm matched Tracks and Clusters
        Map<Track, Cluster> matchedTrackClusterMap = new HashMap<Track,Cluster>();
        //Current (04.19.21) matcher requires List of List of Tracks
        List<List<Track>> trackCollections = new ArrayList<List<Track>>();
        trackCollections.add(tracks);
        //Run Track Cluster Matcher
        matchedTrackClusterMap = matcher.matchTracksToClusters(event, trackCollections, clusters, cuts, -1, false, true, ecal, beamEnergy);
        //Initialize variables for counting instances
        int nele = 0;
        int npos = 0;
        int nphoton = 0;
        int nGoodMatches = 0;
        int nBadMatches = 0;
        int nMissedMatches = 0;
        int nUnmatchedAmbiguousClusterInefficiency = 0;
        int nMatchedAmbiguousTrackInefficiency = 0;
        int nMatchedAmbiguousClusterInefficiency = 0;
        int nUnmatchedAmbiguousTrackInefficiency = 0;
        int nPhotonsMatchedToTrack = 0;
        int nTruthPhotons = 0;
        int nRogueMatches = 0;
        int nIddPhotons = 0;
        int nRogueTracks = 0;
        int nTracklessMCPClusterMatchedToRogueTrack = 0;
        int nTracklessMCPCorrectlyUnmatched = 0;
        int nClusterlessMCPMatchedToRogueCluster = 0;
        int nClusterlessMCPCorrectlyUnmatched = 0;
        int nMCParticlesEvaluated = 0;

        //Initialize Maps/Lists for binning each case
        //Each case below is orthogonal, and all entries will add up to the
        //total number of MCPs in the mcpFullTruthMap
        Map<MCParticle,Pair<Track,Cluster>> goodMatches = new HashMap<MCParticle,Pair<Track,Cluster>>();
        Map<MCParticle,Pair<Track,Cluster>> badMatches = new HashMap<MCParticle,Pair<Track,Cluster>>();
        Map<MCParticle,Pair<Track,Cluster>> unmatchedAmbiguousClusterInefficiency = new HashMap<MCParticle,Pair<Track,Cluster>>();
        Map<MCParticle,Pair<Track,Cluster>> photonsMatchedToTrack = new HashMap<MCParticle,Pair<Track,Cluster>>();
        Map<MCParticle,Pair<Track,Cluster>> truthPhotonClusters = new HashMap<MCParticle,Pair<Track,Cluster>>();
        Map<MCParticle,Pair<Track,Cluster>> matchedAmbiguousTrackInefficiency = new HashMap<MCParticle,Pair<Track,Cluster>>();
        Map<MCParticle,Pair<Track,Cluster>> matchedAmbiguousClusterInefficiency = new HashMap<MCParticle,Pair<Track,Cluster>>();
        Map<MCParticle, Pair<Track,Cluster>> unmatchedAmbiguousTrackInefficiency = new HashMap<MCParticle, Pair<Track,Cluster>>();
        Map<Track, Cluster> rogueMatches = new HashMap<Track,Cluster>();
        List<Track> rogueTracks = new ArrayList<Track>();
        List<Cluster> iddPhotons = new ArrayList<Cluster>();
        List<MCParticle> mcparticlesEvaluated = new ArrayList<MCParticle>();

        //EVALUTE PERFORMANCE OF MATCHER USING TRUTH INFORMATION 
        //Loop over the mcpFullTruthMap and compare truth Track+Cluster pairing
        //with the Track+Cluster pairs returned by the matcher.
        for(Map.Entry<MCParticle, Pair<Cluster, List<Track>>> entry : mcpFullTruthMap.entrySet()){
            MCParticle mcp = entry.getKey();
            int nmcpHits = getNSimTrackerHits(event, mcp);
            Pair<Cluster, List<Track>> pair = entry.getValue();
            Cluster mcpCluster = pair.getFirstElement();
            List<Track> mcpTracks = pair.getSecondElement();

            //Check charged MCP's first
            if(Math.abs(mcp.getPDGID()) == 11){

                if(mcp.getPDGID() == 11){
                    nele = nele + 1;
                }
                if(mcp.getPDGID() == -11){
                    npos = npos + 1;
                }

                //If MCP has truthCluster
                //Check matcher algorithm results to see if truthCluster is matched to
                //algorithmTrack.
                if(mcpCluster != null){

                    //If truthCluster is matched to algorithmTrack
                    if(matchedTrackClusterMap.containsValue(mcpCluster)){
                        //get algorithmTrack
                        Track algTrack = null;
                        for(Map.Entry<Track, Cluster> subentry : matchedTrackClusterMap.entrySet()){
                            if(subentry.getValue() == mcpCluster){
                                algTrack = subentry.getKey();
                                break;
                            }
                        }

                        double[] trackP = TrackUtils.getTrackStateAtLocation(algTrack,TrackState.AtIP).getMomentum();
                        double trackPmag = Math.sqrt(Math.pow(trackP[0],2) + Math.pow(trackP[1],2) + Math.pow(trackP[2],2));

                        //If algorithmTrack has truth information...
                        //Does algorithmTrackMCP = truthClusterMCP?
                        if(trackMCPMapDisamb.containsKey(algTrack)){
                            //If yes, Track+Cluster matching algorithm has
                            //correctly matched a Track and Cluster together,
                            //based on available truth information. 
                            if(mcpTracks.contains(algTrack)){
                                goodMatches.put(mcp, new Pair<Track, Cluster>(algTrack, mcpCluster));
                                nGoodMatches = nGoodMatches + 1;
                            }
                            //If no, matching algorithm has matched Cluster to
                            //the wrong Track, according to truth information.
                            else{
                                badMatches.put(mcp, new Pair<Track, Cluster>(algTrack, mcpCluster));
                                nBadMatches = nBadMatches + 1;
                                plots2D.get("mcp_loop_truth_cluster_matched_to_wrong_track_p_v_p").fill(mcp.getMomentum().magnitude(),trackPmag);
                            }
                        }

                        //If algorithmTrack has no available truth information
                        //We can't necessarily validate this truthCluster +
                        //algTrack match returned by the matcher
                        else{
                            //If algTrack has no truth information, and
                            //was matched to cluster of a MCP that does
                            //have truth tracks, then we know this algTrack
                            //match is wrong.
                            if(mcpTracks.size() > 0){
                                badMatches.put(mcp, new Pair<Track, Cluster>(algTrack, mcpCluster));
                                nBadMatches = nBadMatches + 1;
                                plots2D.get("mcp_loop_truth_cluster_matched_to_wrong_track_p_v_p").fill(mcp.getMomentum().magnitude(),trackPmag);
                            }
                            //If MCP has no truth track, we cant be sure that
                            //this algTrack is wrong
                            else{
                                //The absence of truthTracks indicates either
                                //1. Tracking Inefficiency
                                //2. Track Truth Matching Inefficiency 
                                //3. MCP is outside Tracker acceptance, and
                                //leaves no reconstructable Track
                                
                                //If MCP leaves < 6 SimTrackerHits, MCP does
                                //not leave a reconstructable Track
                                if(getNSimTrackerHits(event,mcp) < 6){
                                    nTracklessMCPClusterMatchedToRogueTrack = nTracklessMCPClusterMatchedToRogueTrack + 1;
                                }
                                //Else there is some ambiguous tracking
                                //inefficiency 
                                else{
                                    matchedAmbiguousTrackInefficiency.put(mcp, new Pair<Track, Cluster>(algTrack, mcpCluster));
                                    nMatchedAmbiguousTrackInefficiency = nMatchedAmbiguousTrackInefficiency + 1; 
                                    plots1D.get("mcp_loop_trackless_truth_cluster_matched_to_rogue_track_nSimTrackerHits").fill(nmcpHits);
                                }
                            }
                        }
                    }

                    //If mcpCluster is not matched by the Matcher to any Track
                    else{

                        //If MCP has mcpTracks to match to, then the Matching
                        //algorithm has missed this match
                        if(mcpTracks.size() > 0){
                            //missedMatches.put(mcp, new Pair<Track, Cluster>(null, mcpCluster));
                            nMissedMatches = nMissedMatches + 1;
                        }

                        //If MCP has no mcpTracks...
                        else{
                            //If MCP does not leave a reconstructable
                            //Track...it's good that the Matching algorithm did
                            //not match its mcpCluster to any Track
                            if(getNSimTrackerHits(event,mcp) < 6){
                                //tracklessMCPClusterNotMatchedToTrack.put(mcp, new Pair<null, mcpCluster>);
                                nTracklessMCPCorrectlyUnmatched = nTracklessMCPCorrectlyUnmatched + 1;
                            }

                            //If MCP leaves > 5 SimTrackerHits, its truthTrack
                            //should have been found...but wasn't because of
                            //some inefficiency
                            else{
                                unmatchedAmbiguousTrackInefficiency.put(mcp, new Pair<Track, Cluster>(null, mcpCluster));
                                nUnmatchedAmbiguousTrackInefficiency = nUnmatchedAmbiguousTrackInefficiency + 1;
                                plots1D.get("mcp_loop_no_truth_tracks_truth_cluster_not_matched_to_track_nSimTrackerHits").fill(nmcpHits);
                            }
                        }
                    }
                }

                //If MCP has NO truthCluster...
                //We check if any of the mcpTracks were matched to any Clusters
                else{
                    //If MCP leaves 0 SimCalHits, we know that there is no
                    //truthCluster for this MCP, so we don't expect this
                    //mcpTrack to be matached to any Clusters

                    boolean clusterExpected = false;
                    boolean algClusterIsRogue = true;


                    //loop over all mcpTracks on this MCP
                    Cluster matchedCluster = null;
                    Track matchedTrack = null;
                    //get mcpTrack with most hits
                    Track mcpBestTrack = getMcpBestTrack(mcp, mcpTrackMapDisamb);
                    //Check if mcpTrack was matched to Cluster by Matcher
                    matchedCluster = matchedTrackClusterMap.get(mcpBestTrack);
                    matchedTrack = mcpBestTrack;
                    //If best mcpTrack wasn't matched to a Cluster, see if any
                    //of the lesser mcpTracks were
                    if(matchedCluster == null){
                        for(Track track : mcpTracks){
                            matchedCluster = matchedTrackClusterMap.get(track);
                            if(matchedCluster != null){
                                matchedTrack = track;
                                break;
                            }
                        }
                    }

                    //If MCP > 0 SimCalHits, *maybe* there should be a
                    //mcpCluster.
                    if(getNSimCalHits(event, "EcalHits", mcp) > 0){
                        clusterExpected = true;
                    }

                    //If MCP possibly left a Cluster that wasn't found using
                    //truth info...
                    if(clusterExpected){
                        //If mcpTrack unmatched to Cluster, this case is
                        //ambiguous, as the truth of whether or not this MCP
                        //left a Cluster is unknown
                        if(matchedCluster == null){
                            unmatchedAmbiguousClusterInefficiency.put(mcp, new Pair<Track, Cluster>(matchedTrack, null));
                            nUnmatchedAmbiguousClusterInefficiency = nUnmatchedAmbiguousClusterInefficiency + 1;
                        }

                        //if matchedCluster belongs to some other MCP, we know
                        //this match is wrong
                        else if(clustersMCPMap.containsKey(matchedCluster)){
                            badMatches.put(mcp, new Pair<Track, Cluster>(matchedTrack,matchedCluster));
                            nBadMatches = nBadMatches + 1;
                            plots2D.get("mcp_loop_truth_track_matched_to_wrong_cluster_E_v_E").fill(mcp.getEnergy(),matchedCluster.getEnergy());
                        }

                        //If matchedCluster has no truth info, we can't say if
                        //this mcpTrack + matchedCluster pair is wrong or not
                        else{
                            matchedAmbiguousClusterInefficiency.put(mcp, new Pair<Track, Cluster>(matchedTrack, matchedCluster));
                            nMatchedAmbiguousClusterInefficiency = nMatchedAmbiguousClusterInefficiency + 1; 
                        }
                    }

                    //If MCP leaves no Cluster...
                    else{

                        //mcpTrack should not be matched to any Cluster
                        if(matchedCluster == null){
                            //clusterlessMCPtrackNotMatchedToCluster.put(mcp, new Pair<Track, Cluster>(matchedTrack, null));
                            nClusterlessMCPCorrectlyUnmatched = nClusterlessMCPCorrectlyUnmatched + 1;
                        }

                        //If mcpTrack matched to otherMcpCluster, bad match
                        else if(clustersMCPMap.containsKey(matchedCluster)){
                            badMatches.put(mcp, new Pair<Track, Cluster>(matchedTrack,matchedCluster));
                            nBadMatches = nBadMatches + 1;
                            plots2D.get("mcp_loop_truth_track_matched_to_wrong_cluster_E_v_E").fill(mcp.getEnergy(),matchedCluster.getEnergy());
                        }

                        //If clusterless mcpTrack matched to truthless Cluster
                        else{
                            nClusterlessMCPMatchedToRogueCluster = nClusterlessMCPMatchedToRogueCluster + 1;
                        }
                    }
                }
            }

            //Loop over photon MCPs
            if(Math.abs(mcp.getPDGID()) == 22){
                nphoton = nphoton + 1;
                //If photon is matched to any track            
                Track matchedTrack = null;
                if(matchedTrackClusterMap.containsValue(mcpCluster)){
                    for(Map.Entry<Track, Cluster> subentry : matchedTrackClusterMap.entrySet()){
                        if(subentry.getValue() == mcpCluster){
                            matchedTrack = subentry.getKey();
                            break;
                        }
                    }
                    photonsMatchedToTrack.put(mcp, new Pair<Track, Cluster>(matchedTrack, mcpCluster));
                    nPhotonsMatchedToTrack = nPhotonsMatchedToTrack + 1;
                }
                else{
                    truthPhotonClusters.put(mcp, new Pair<Track, Cluster>(null, mcpCluster));
                    nTruthPhotons = nTruthPhotons + 1;
                }
            }
            mcparticlesEvaluated.add(mcp);
            nMCParticlesEvaluated = nMCParticlesEvaluated + 1;
        }

        plots1D.get("mcp_loop_nEle_per_event").fill(nele);
        plots1D.get("mcp_loop_nPos_per_event").fill(npos);
        plots1D.get("mcp_loop_nPhotons_per_event").fill(nphoton);
        plots1D.get("mcp_loop_nCharged_particles_per_event").fill(nele + npos);
        plots2D.get("mcp_loop_nCharged_particles_v_nTracks_per_event").fill(nele + npos, tracks.size());
        plots2D.get("mcp_loop_nMCPs_v_nClusters_per_event").fill(nele + npos + nphoton,clusters.size());


        //Loop over rogue clusters to find rogue matches and rogue photons
        for(Cluster cluster : clusters){
            if(clustersMCPMap.containsKey(cluster))
                continue;
            if(matchedTrackClusterMap.containsValue(cluster)){
                Track matchedTrack = null;
                for(Map.Entry<Track, Cluster> subentry : matchedTrackClusterMap.entrySet()){
                    if(subentry.getValue() == cluster){
                        matchedTrack = subentry.getKey();
                    }
                }
                if(trackMCPMapDisamb.containsKey(matchedTrack)){
                    continue;
                }
                else{
                    rogueMatches.put(matchedTrack, cluster);
                    nRogueMatches = nRogueMatches + 1;
                }
            }
            else{
                iddPhotons.add(cluster);
                nIddPhotons = nIddPhotons + 1;
            }
        }

        //Loop over rogue tracks to find rogue matches and rogue tracks
        for(Track track : tracks){
            if(trackMCPMapDisamb.containsKey(track))
                continue;
            Cluster matchedCluster = matchedTrackClusterMap.get(track);
            if(matchedCluster == null){
                rogueTracks.add(track);
                nRogueTracks = nRogueTracks + 1;
            }
            if(clustersMCPMap.containsKey(matchedCluster))
                continue;
            else
                if(!rogueMatches.containsKey(track)){
                    rogueMatches.put(track, matchedCluster);
                    nRogueMatches = nRogueMatches + 1;
                }
        }

        //Fill count histograms
        plots1D.get("nGoodMatches").fill(nGoodMatches);
        plots1D.get("nGoodMatches").fill(-1, nGoodMatches);
        goodMatchPlots(goodMatches);

        plots1D.get("nBadMatches").fill(nBadMatches);
        plots1D.get("nBadMatches").fill(-1, nBadMatches);

        plots1D.get("nMissedMatches").fill(nMissedMatches);
        plots1D.get("nMissedMatches").fill(-1, nMissedMatches);

        plots1D.get("nUnmatchedAmbiguousClusterInefficiency").fill(nUnmatchedAmbiguousClusterInefficiency);
        plots1D.get("nUnmatchedAmbiguousClusterInefficiency").fill(-1, nUnmatchedAmbiguousClusterInefficiency);

        plots1D.get("nMatchedAmbiguousTrackInefficiency").fill(nMatchedAmbiguousTrackInefficiency);
        plots1D.get("nMatchedAmbiguousTrackInefficiency").fill(-1, nMatchedAmbiguousTrackInefficiency);

        plots1D.get("nMatchedAmbiguousClusterInefficiency").fill(nMatchedAmbiguousClusterInefficiency);
        plots1D.get("nMatchedAmbiguousClusterInefficiency").fill(-1, nMatchedAmbiguousClusterInefficiency);

        plots1D.get("nUnmatchedAmbiguousTrackInefficiency").fill(nUnmatchedAmbiguousTrackInefficiency);
        plots1D.get("nUnmatchedAmbiguousTrackInefficiency").fill(-1, nUnmatchedAmbiguousTrackInefficiency);

        plots1D.get("nPhotons_matched_to_track").fill(nPhotonsMatchedToTrack);
        plots1D.get("nPhotons_matched_to_track").fill(-1, nPhotonsMatchedToTrack);

        plots1D.get("nTruth_photon_clusters").fill(nTruthPhotons);
        plots1D.get("nTruth_photon_clusters").fill(-1, nTruthPhotons);

        plots1D.get("nRogueMatches").fill(nRogueMatches);
        plots1D.get("nRogueMatches").fill(-1, nRogueMatches);

        plots1D.get("nRogueTracks").fill(nRogueTracks);
        plots1D.get("nRogueTracks").fill(-1, nRogueTracks);

        plots1D.get("nIddPhotons").fill(nIddPhotons);
        plots1D.get("nIddPhotons").fill(-1, nIddPhotons);

        plots1D.get("nTracklessMCPClusterMatchedToRogueTrack").fill(nTracklessMCPClusterMatchedToRogueTrack);
        plots1D.get("nTracklessMCPClusterMatchedToRogueTrack").fill(-1, nTracklessMCPClusterMatchedToRogueTrack);

        plots1D.get("nTracklessMCPCorrectlyUnmatched").fill(nTracklessMCPCorrectlyUnmatched);
        plots1D.get("nTracklessMCPCorrectlyUnmatched").fill(-1, nTracklessMCPCorrectlyUnmatched);

        plots1D.get("nClusterlessMCPMatchedToRogueCluster").fill(nClusterlessMCPMatchedToRogueCluster);
        plots1D.get("nClusterlessMCPMatchedToRogueCluster").fill(-1, nClusterlessMCPMatchedToRogueCluster);

        plots1D.get("nClusterlessMCPCorrectlyUnmatched").fill(nClusterlessMCPCorrectlyUnmatched);
        plots1D.get("nClusterlessMCPCorrectlyUnmatched").fill(-1, nClusterlessMCPCorrectlyUnmatched);

        plots1D.get("nMCParticlesEvaluated").fill(nMCParticlesEvaluated);
        plots1D.get("nMCParticlesEvaluated").fill(-1, nMCParticlesEvaluated);
        */
    }

    private void drawEcalFace(List<Cluster> clusters){
        //Define line that draws the beam gap Ecal crystal edge
        int nx = 46;
        int ny = 5;
        double crystalface = 13.0; //mm
        
        for(Cluster cluster : clusters){

            double clusterx = cluster.getPosition()[0];
            double clustery = cluster.getPosition()[1];
            double clusterz = cluster.getPosition()[2];
            
            double leftx = clusterx - (crystalface/2);
            double rightx = clusterx + (crystalface/2);
            double upy = clustery + (crystalface/2);
            double downy = clustery - (crystalface/2);

            plots2D.get("ecal_crystal_acceptance_xy").fill(clusterx,clustery);
            plots2D.get("ecal_crystal_acceptance_xy").fill(leftx,clustery);
            plots2D.get("ecal_crystal_acceptance_xy").fill(rightx,clustery);
            plots2D.get("ecal_crystal_acceptance_xy").fill(clusterx,upy);
            plots2D.get("ecal_crystal_acceptance_xy").fill(clusterx,downy);

            plots2D.get("ecal_crystal_acceptance_xy").fill(rightx,upy);
            plots2D.get("ecal_crystal_acceptance_xy").fill(rightx,downy);
            plots2D.get("ecal_crystal_acceptance_xy").fill(leftx,downy);
            plots2D.get("ecal_crystal_acceptance_xy").fill(leftx,upy);
        }
    }

    private void disambiguateTruthTracks(Map<MCParticle, Map<Track, Integer>> mcpTrackMap, Map<MCParticle, Map<Track, Integer>> mcpTrackMapDisamb , Map<Track, Map<MCParticle, Integer>> trackMCPMap, Map<Track, MCParticle> trackMCPMapDisamb){ 

        for(Map.Entry<Track, Map<MCParticle, Integer>> entry : trackMCPMap.entrySet()){

            if(entry.getValue() == null){
                continue;
            }

            Track track = entry.getKey();
            double trackPmag = new BasicHep3Vector(track.getTrackStates().get(0).getMomentum()).magnitude();
            int mosthits = 0;
            MCParticle bestMCP = null;
            double bestMCPmomentum = 0.0;
            Map<MCParticle, Integer> bMap = entry.getValue();

            for(Map.Entry<MCParticle, Integer> subentry : bMap.entrySet()){
                //If Track is tied with two MCPs, match to Track with momentum
                //nearest MCP
                if(subentry.getValue() == mosthits){
                    if(Math.abs((1-(subentry.getKey().getMomentum().magnitude()/trackPmag))) <= 
                            Math.abs((1-(bestMCPmomentum/trackPmag)))){
                        bestMCP = subentry.getKey();
                        bestMCPmomentum = subentry.getKey().getMomentum().magnitude();
                        mosthits = subentry.getValue();
                    }
                }
                //Match Track to MCP that leaves most hits on Track
                if(subentry.getValue() > mosthits){
                    mosthits = subentry.getValue();
                    bestMCP = subentry.getKey();
                    bestMCPmomentum = subentry.getKey().getMomentum().magnitude();
                }
            }

            //Add track and best MCP to map
            //trackMCPMapDisamb.put(entry.getKey(),bestMCP);
            trackMCPMapDisamb.put(track, bestMCP);

            //Fill disambiguated map
            if(!mcpTrackMapDisamb.containsKey(bestMCP)){
                Map<Track, Integer> tmpMap = new HashMap<Track,Integer>();
                tmpMap.put(entry.getKey(),mosthits);
                mcpTrackMapDisamb.put(bestMCP,tmpMap);
            }
            else{
                Map<Track, Integer> tmpMap = mcpTrackMapDisamb.get(bestMCP);
                tmpMap.put(entry.getKey(),mosthits);
                mcpTrackMapDisamb.put(bestMCP,tmpMap);
            }
        }
    }

    private void buildTrackTruthMaps(List<Track> tracks, RelationalTable rawtomc, Map<MCParticle, 
            Map<Track, Integer>> mcpTrackMapDisamb, Map<Track, MCParticle> trackMCPMapDisamb){

        //Track truth matching
        Map<MCParticle, Map<Track, Integer>> mcpTrackMap = new HashMap<MCParticle, Map<Track, Integer>>();
        Map<Track, Map<MCParticle, Integer>> trackMCPMap = new HashMap<Track, Map<MCParticle, Integer>>();
        for(Track track : tracks){
            double trackPmag = new BasicHep3Vector(track.getTrackStates().get(0).getMomentum()).magnitude();
            int charge = -1* (int)Math.signum(track.getTrackStates().get(0).getOmega());

            //run Track to MCP truth matching
            TrackTruthMatcher truthMatcher = new TrackTruthMatcher(track, rawtomc, 0, 6);
            
            //Get all MCPs that leave hits on this Track
            Map<MCParticle, Set<Integer>> mcpsOnTrack = truthMatcher.getLayerHitsForAllMCPs();

            if(mcpsOnTrack == null)
                continue;

            for(Map.Entry<MCParticle, Set<Integer>> entry : mcpsOnTrack.entrySet()){
                MCParticle mcp = entry.getKey();
                int mcp_charge = (int) mcp.getCharge();
                int nhits = entry.getValue().size();

                //Only match a track to a MCP if it has > 6 MCP hits
                if(nhits < 6)
                    continue;

                plots2D.get("track_matched_to_mcps_p_v_p").fill(trackPmag*charge,mcp.getMomentum().magnitude()*mcp_charge); 

                //Fill Track to MCP Map
                if(!trackMCPMap.containsKey(track)){
                    Map<MCParticle, Integer> mcpMap = new HashMap<MCParticle, Integer>();
                    mcpMap.put(mcp, nhits);
                    trackMCPMap.put(track, mcpMap);
                }
                else{
                    Map<MCParticle, Integer> mcpMap = trackMCPMap.get(track);
                    mcpMap.put(mcp, nhits);
                    trackMCPMap.put(track, mcpMap);
                }

                //Fill MCP to Track Map
                if(!mcpTrackMap.containsKey(mcp)){
                    Map<Track, Integer> trackMap = new HashMap<Track, Integer>();
                    trackMap.put(track, nhits);
                    mcpTrackMap.put(mcp, trackMap);
                }

                else{
                    Map<Track, Integer> trackMap = mcpTrackMap.get(mcp);
                    trackMap.put(track, nhits);
                    mcpTrackMap.put(mcp, trackMap);
                }
            }

            MCParticle bestmcp = truthMatcher.getMCParticle();
            int mcp_charge = (int) bestmcp.getCharge();
            int bestNhits = truthMatcher.getNGoodHits();
            if(truthMatcher.getMCParticle() != null){
                plots2D.get("track_matched_to_best_mcp_p_v_p").fill(trackPmag*charge,bestmcp.getMomentum().magnitude()*mcp_charge); 
                plots1D.get("track_matched_to_best_mcp_nHits").fill(bestNhits); 
            }

            plots1D.get("track_matched_to_n_mcps").fill(trackMCPMap.get(track).size());
        }

        System.out.println("Number of MCPs matched to track new: " + mcpTrackMap.size());
        //Ensure that a track is only matched to a unique MCParticle by
        //disambiguating matches
        //Map<MCParticle, Map<Track, Integer>> mcpTrackMapDisamb = new HashMap<MCParticle, Map<Track, Integer>>();
        //Map<Track, MCParticle> trackMCPMapDisamb = new HashMap<Track, MCParticle>();
        disambiguateTruthTracks(mcpTrackMap, mcpTrackMapDisamb, trackMCPMap, trackMCPMapDisamb);
        System.out.println("Number of MCPs matched to track new after disamb" + mcpTrackMapDisamb.size());
    }


    private void trackClusterResidualParameterization(Map<MCParticle,Pair<Track, Cluster>> mcpTrackClusterPairs){

        //baboon
        for(Map.Entry<MCParticle, Pair<Track, Cluster>> entry : mcpTrackClusterPairs.entrySet()){

            Track track = entry.getValue().getFirstElement();
            Cluster cluster = entry.getValue().getSecondElement();
            int charge = -1* (int) Math.signum(track.getTrackStates().get(0).getOmega());
            double[] trackP = TrackUtils.getTrackStateAtLocation(track,TrackState.AtIP).getMomentum();
            double trackPmag = new BasicHep3Vector(track.getTrackStates().get(0).getMomentum()).magnitude();
            double [] params = track.getTrackParameters();
            double tanlambda = params[4];
            boolean isTop;

            List<Double> trackPos = getTrackPositionAtEcal(track);
            double trackx = trackPos.get(0);
            double tracky = trackPos.get(1);
            double trackz = trackPos.get(2);

            double clusterx = cluster.getPosition()[0];
            double clustery = cluster.getPosition()[1];
            double clusterz = cluster.getPosition()[2];
            double clusterEnergy = cluster.getEnergy();

            if(charge < 0){
                if(tanlambda > 0){
                    plots2D.get(String.format("%s_ele_TOP_track_cluster_truth_pairs_residual_param_dx",this.trackCollectionName)).fill(trackPmag,clusterx-trackx);
                    plots2D.get(String.format("%s_ele_TOP_track_cluster_truth_pairs_residual_param_dy",this.trackCollectionName)).fill(trackPmag,clustery-tracky);
                    plots2D.get(String.format("%s_ele_TOP_track_cluster_truth_pairs_residual_param_dz",this.trackCollectionName)).fill(trackPmag,clusterz-trackz);
                }
                else{
                    plots2D.get(String.format("%s_ele_BOTTOM_track_cluster_truth_pairs_residual_param_dx",this.trackCollectionName)).fill(trackPmag,clusterx-trackx);
                    plots2D.get(String.format("%s_ele_BOTTOM_track_cluster_truth_pairs_residual_param_dy",this.trackCollectionName)).fill(trackPmag,clustery-tracky);
                    plots2D.get(String.format("%s_ele_BOTTOM_track_cluster_truth_pairs_residual_param_dz",this.trackCollectionName)).fill(trackPmag,clusterz-trackz);
                }
            }
            else{
                if(tanlambda > 0){
                    plots2D.get(String.format("%s_pos_TOP_track_cluster_truth_pairs_residual_param_dx",this.trackCollectionName)).fill(trackPmag,clusterx-trackx);
                    plots2D.get(String.format("%s_pos_TOP_track_cluster_truth_pairs_residual_param_dy",this.trackCollectionName)).fill(trackPmag,clustery-tracky);
                    plots2D.get(String.format("%s_pos_TOP_track_cluster_truth_pairs_residual_param_dz",this.trackCollectionName)).fill(trackPmag,clusterz-trackz);
                }
                else{
                    plots2D.get(String.format("%s_pos_BOTTOM_track_cluster_truth_pairs_residual_param_dx",this.trackCollectionName)).fill(trackPmag,clusterx-trackx);
                    plots2D.get(String.format("%s_pos_BOTTOM_track_cluster_truth_pairs_residual_param_dy",this.trackCollectionName)).fill(trackPmag,clustery-tracky);
                    plots2D.get(String.format("%s_pos_BOTTOM_track_cluster_truth_pairs_residual_param_dz",this.trackCollectionName)).fill(trackPmag,clusterz-trackz);
                }

            }

        }

    }

    private void goodMatchPlots(Map<MCParticle,Pair<Track, Cluster>> goodMatches){

        for(Map.Entry<MCParticle, Pair<Track, Cluster>> entry : goodMatches.entrySet()){

            Track track = entry.getValue().getFirstElement();
            Cluster cluster = entry.getValue().getSecondElement();
            int charge = -1* (int) Math.signum(track.getTrackStates().get(0).getOmega());
            double[] trackP = TrackUtils.getTrackStateAtLocation(track,TrackState.AtIP).getMomentum();
            double trackPmag = new BasicHep3Vector(track.getTrackStates().get(0).getMomentum()).magnitude();
            double [] params = track.getTrackParameters();
            double tanlambda = params[4];
            boolean isTop;

            List<Double> trackPos = getTrackPositionAtEcal(track);
            double trackx = trackPos.get(0);
            double tracky = trackPos.get(1);
            double trackz = trackPos.get(2);

            double clusterx = cluster.getPosition()[0];
            double clustery = cluster.getPosition()[1];
            double clusterz = cluster.getPosition()[2];
            double clusterEnergy = cluster.getEnergy();

            if(charge < 0){
                plots2D.get(String.format("%s_ele_good_matches_trackP_v_clusterE",this.trackCollectionName)).fill(trackPmag, cluster.getEnergy());
                if(tanlambda > 0){
                    plots2D.get(String.format("%s_ele_TOP_track_cluster_good_matches_dx",this.trackCollectionName)).fill(trackPmag,clusterx-trackx);
                    plots2D.get(String.format("%s_ele_TOP_track_cluster_good_matches_dy",this.trackCollectionName)).fill(trackPmag,clustery-tracky);
                }
                else{
                    plots2D.get(String.format("%s_ele_BOTTOM_track_cluster_good_matches_dx",this.trackCollectionName)).fill(trackPmag,clusterx-trackx);
                    plots2D.get(String.format("%s_ele_BOTTOM_track_cluster_good_matches_dy",this.trackCollectionName)).fill(trackPmag,clustery-tracky);
                }
            }
            else{
                plots2D.get(String.format("%s_pos_good_matches_trackP_v_clusterE",this.trackCollectionName)).fill(trackPmag, cluster.getEnergy());
                if(tanlambda > 0){
                    plots2D.get(String.format("%s_pos_TOP_track_cluster_good_matches_dx",this.trackCollectionName)).fill(trackPmag,clusterx-trackx);
                    plots2D.get(String.format("%s_pos_TOP_track_cluster_good_matches_dy",this.trackCollectionName)).fill(trackPmag,clustery-tracky);
                }
                else{
                    plots2D.get(String.format("%s_pos_BOTTOM_track_cluster_good_matches_dx",this.trackCollectionName)).fill(trackPmag,clusterx-trackx);
                    plots2D.get(String.format("%s_pos_BOTTOM_track_cluster_good_matches_dy",this.trackCollectionName)).fill(trackPmag,clustery-tracky);
                }
            }
        }
    }

    public boolean isTrackInEcal(double trackx, double tracky){

        //Define first order Ecal geometry --> assumed square here,
        //smaller than actual beamgap. Improve x geometry 
        double ecalx1 = -276.0; //mm
        double ecalx2 = 361.0;
        double ecaly1 = 91.0;
        double ecaly2 = -91.0;
        double bgapup = 22;
        double bgapdown = -22;

        double eholex11 = -93.0;
        double eholex12 = -70.0;

        double eholex22 = 15.0;
        double eholex21 = 29;

        double eholey12 = 36.0;
        double eholey11 = 22.4;

        double eholey22 = -36.0;
        double eholey21 = -22.3;

        boolean inEcalAccept = true;

        if(trackx < ecalx1 || trackx > ecalx2)
            inEcalAccept = false;
        if(tracky > ecaly1 || tracky < ecaly2)
            inEcalAccept = false;
        if(tracky < bgapup && tracky > bgapdown)
            inEcalAccept = false;
        if((trackx > eholex12 && trackx < eholex22) && ( (tracky < eholey12) && (tracky > eholey11)))
            inEcalAccept = false;
        if((trackx > eholex12 && trackx < eholex22) && ( (tracky < eholey21) && (tracky > eholey22)))
            inEcalAccept = false;
        return inEcalAccept;
    }


    public SimTrackerHit getTrackScoringPlaneHit(EventHeader event, MCParticle mcp, String ecalScoringPlaneHitsCollectionName) {

        
        //If even doesnt have collection of Ecal scoring plane hits, skip
        if(!event.hasCollection(SimTrackerHit.class, ecalScoringPlaneHitsCollectionName)) return null;

        List<SimTrackerHit> scoringPlaneHits = event.get(SimTrackerHit.class, ecalScoringPlaneHitsCollectionName);


        //Check for simtrackerhit MCP that matches trackMCP
        if(mcp == null)
            return null;
        SimTrackerHit matchedScoringPlaneHit = null;
        for(SimTrackerHit scoringPlaneHit : scoringPlaneHits){
            // If the MC particles don't match, move on to the next particle
            if(!(scoringPlaneHit.getMCParticle() == mcp)) continue;
            matchedScoringPlaneHit = scoringPlaneHit;
            // Once a match is found, there is no need to loop through the rest of the list
            break;
        }
        return matchedScoringPlaneHit;
    }
    

    //Has bug, doesnt work
    public double[] getExtrapolatedTrackScoringPlaneHit(EventHeader event, Track track, SimTrackerHit scoringplaneHit){
    
        double truthxpos;
        double truthypos;
        double truthzpos;

        double trackT;
        double simTrackT;

        simTrackT = scoringplaneHit.getTime();

        truthxpos = scoringplaneHit.getPoint()[0];
        truthypos = scoringplaneHit.getPoint()[1];
        truthzpos = scoringplaneHit.getPoint()[2];

        //multiply charge by factor of -1 (WHY!?)
        int charge = -1* (int) Math.signum(track.getTrackStates().get(0).getOmega());

        double[] truthP = scoringplaneHit.getMomentum();
        double truthPmag = Math.sqrt(Math.pow(truthP[0],2) + Math.pow(truthP[1],2) + Math.pow(truthP[2],2));

        //In PhysicsRun2016 detector geometry, scoring plane is located
        //upstream of Ecal face by ~30 mm...
        //Need to extrapolate MCP scoring plane hits to the Ecal, so that we
        //can compare track positions at Ecal with MCP hit at Ecal
        //rotate HPS coordinates to Tracking coordinates (XYZ)HPS -> (ZXY)
        Hep3Vector truthPVecTracking = new BasicHep3Vector(truthP[2],truthP[0],truthP[1]);
        Hep3Vector truthPosVecTracking = new BasicHep3Vector(truthzpos,truthxpos,truthypos);

        Pair<Hep3Vector,Hep3Vector> extrapVecsTracking = null;
        Detector detector = event.getDetector();
        fM = detector.getFieldMap();
        RK4integrator RKint = new RK4integrator(charge,1, fM);
        extrapVecsTracking = RKint.integrate(truthPosVecTracking,truthPVecTracking,1393-scoringplaneHit.getPoint()[2]);
        Hep3Vector RKextrapPos = extrapVecsTracking.getFirstElement();
        truthxpos = RKextrapPos.y();
        truthypos = RKextrapPos.z();
        truthzpos = RKextrapPos.x();
        double[] truthpos = {truthxpos, truthypos, truthzpos};

        return truthpos;
    }


    public void trackScoringPlanePlots(EventHeader event, Track track, SimTrackerHit scoringplaneHit) {

        //Comparison plots for Track extrapolated to ECal vs Truth ScoringPlaneHit at Ecal
    
        double trackx;
        double tracky;
        double trackz;
        double truthxpos;
        double truthypos;
        double truthzpos;
        double dxoffset;

        double trackT;
        double simTrackT;

        if (this.trackCollectionName.contains("GBLTracks")){
            trackT = TrackUtils.getTrackTime(track, hitToStrips, hitToRotated);
        }
        else {
            TrackData trackdata = (TrackData) TrktoData.from(track);
            trackT = trackdata.getTrackTime();
        }
        simTrackT = scoringplaneHit.getTime();

        //Make histograms of truth vs extrapolation
        if(this.trackCollectionName.contains("GBLTracks")) {
            trackx = TrackUtils.getTrackStateAtECal(track).getReferencePoint()[1];
            tracky = TrackUtils.getTrackStateAtECal(track).getReferencePoint()[2];
            trackz = TrackUtils.getTrackStateAtECal(track).getReferencePoint()[0];
            //dxoffset = -4.1;
            dxoffset = 0.0;
        }

        else {
            TrackState ts_ecal = track.getTrackStates().get(track.getTrackStates().size()-1);
            double[] ts_ecalPos = ts_ecal.getReferencePoint();
            trackx = ts_ecalPos[1];
            tracky = ts_ecalPos[2];
            trackz = ts_ecalPos[0];
            dxoffset = 0.0;
        }

        truthxpos = scoringplaneHit.getPoint()[0];
        truthypos = scoringplaneHit.getPoint()[1];
        truthzpos = scoringplaneHit.getPoint()[2];

        //multiply charge by factor of -1 (WHY!?)
        int charge = -1* (int) Math.signum(track.getTrackStates().get(0).getOmega());

        //double[] trackP = track.getMomentum();
        //momentum is rotated coords (x->z, z->y, y->x)
        double[] trackP = TrackUtils.getTrackStateAtLocation(track,TrackState.AtIP).getMomentum();
        double trackPmag = new BasicHep3Vector(track.getTrackStates().get(0).getMomentum()).magnitude();
        double[] truthP = scoringplaneHit.getMomentum();
        double truthPmag = Math.sqrt(Math.pow(truthP[0],2) + Math.pow(truthP[1],2) + Math.pow(truthP[2],2));


        //In PhysicsRun2016 detector geometry, scoring plane is located
        //upstream of Ecal face by ~30 mm...
        //Need to extrapolate MCP scoring plane hits to the Ecal, so that we
        //can compare track positions at Ecal with MCP hit at Ecal
        //rotate HPS coordinates to Tracking coordinates (XYZ)HPS -> (ZXY)
        Hep3Vector truthPVecTracking = new BasicHep3Vector(truthP[2],truthP[0],truthP[1]);
        Hep3Vector truthPosVecTracking = new BasicHep3Vector(truthzpos,truthxpos,truthypos);

        Pair<Hep3Vector,Hep3Vector> extrapVecsTracking = null;
        Detector detector = event.getDetector();
        fM = detector.getFieldMap();
        RK4integrator RKint = new RK4integrator(charge,1, fM);
        extrapVecsTracking = RKint.integrate(truthPosVecTracking,truthPVecTracking,1393-scoringplaneHit.getPoint()[2]);
        Hep3Vector RKextrapPos = extrapVecsTracking.getFirstElement();
        truthxpos = RKextrapPos.y();
        truthypos = RKextrapPos.z();
        truthzpos = RKextrapPos.x();

        //Take position residuals for track and extrap truth hits at ecal
        double dx = truthxpos - trackx;
        double dy = truthypos - tracky;
        double dz = truthzpos - trackz;
        double dr = Math.sqrt(Math.pow(dx,2) + Math.pow(dy,2) + Math.pow(dz,2));
        double dt = simTrackT - trackT;

        //Make plots
        if(charge < 0) {

            //track momentum. Truth and reco

            //Track X,Y position plots at Ecal
            //Extrapolated Track momentum vs truth position residuals
            plots2D.get(String.format("%s_ele_RK4_scoringplanehit_to_ecal_ZvP",this.trackCollectionName)).fill(truthzpos, trackPmag);
            //Residuals between track extrapolated to Ecal Face, and truth hit
            //extrapolated to Ecal Face
            plots1D.get(String.format("%s_ele_track_scoringplane_hit_dx",this.trackCollectionName)).fill(dx);
            plots1D.get(String.format("%s_ele_track_scoringplane_hit_dy",this.trackCollectionName)).fill(dy);
            plots1D.get(String.format("%s_ele_track_scoringplane_hit_dz",this.trackCollectionName)).fill(dz);
            plots1D.get(String.format("%s_ele_track_scoringplane_hit_dr",this.trackCollectionName)).fill(dr);
            plots1D.get(String.format("%s_ele_track_scoringplane_hit_dt",this.trackCollectionName)).fill(dt);

            plots1D.get(String.format("%s_ele_track_scoringplane_hit_px",this.trackCollectionName)).fill(truthP[0]);
            plots1D.get(String.format("%s_ele_track_scoringplane_hit_py",this.trackCollectionName)).fill(truthP[1]);
            plots1D.get(String.format("%s_ele_track_scoringplane_hit_pz",this.trackCollectionName)).fill(truthP[2]);

            plots2D.get(String.format("%s_ele_scoringplaneHit_v_truth_track_p",this.trackCollectionName)).fill(truthPmag,trackPmag);
        }
        else {

            //Track X,Y position at Ecal
            //Extrapolated Track P vs truth position residuals
            plots2D.get(String.format("%s_pos_RK4_scoringplanehit_to_ecal_ZvP",this.trackCollectionName)).fill(truthzpos, trackPmag);
            //Track vs Cluster residuals at Ecal
            plots1D.get(String.format("%s_pos_track_scoringplane_hit_dx",this.trackCollectionName)).fill(dx);
            plots1D.get(String.format("%s_pos_track_scoringplane_hit_dy",this.trackCollectionName)).fill(dy);
            plots1D.get(String.format("%s_pos_track_scoringplane_hit_dz",this.trackCollectionName)).fill(dz);
            plots1D.get(String.format("%s_pos_track_scoringplane_hit_dr",this.trackCollectionName)).fill(dr);
            plots1D.get(String.format("%s_pos_track_scoringplane_hit_dt",this.trackCollectionName)).fill(dt);

            plots1D.get(String.format("%s_pos_track_scoringplane_hit_px",this.trackCollectionName)).fill(truthP[0]);
            plots1D.get(String.format("%s_pos_track_scoringplane_hit_py",this.trackCollectionName)).fill(truthP[1]);
            plots1D.get(String.format("%s_pos_track_scoringplane_hit_pz",this.trackCollectionName)).fill(truthP[2]);

            plots2D.get(String.format("%s_pos_scoringplaneHit_v_truth_track_p",this.trackCollectionName)).fill(truthPmag,trackPmag);
        }

    }

    public List<Double> getTrackPositionAtEcal(Track track){

        double trackx;
        double tracky;
        double trackz;

        if (this.trackCollectionName.contains("GBLTracks")){
            trackx = TrackUtils.getTrackStateAtECal(track).getReferencePoint()[1];
            tracky = TrackUtils.getTrackStateAtECal(track).getReferencePoint()[2];
            trackz = TrackUtils.getTrackStateAtECal(track).getReferencePoint()[0];
        }
        else {
            TrackData trackdata = (TrackData) TrktoData.from(track);
            TrackState ts_ecal = track.getTrackStates().get(track.getTrackStates().size()-1);
            double[] ts_ecalPos = ts_ecal.getReferencePoint();
            trackx = ts_ecalPos[1];
            tracky = ts_ecalPos[2];
            trackz = ts_ecalPos[0];
        }
        List<Double> trackxyz = new ArrayList<Double>(); 
        trackxyz.add(trackx);
        trackxyz.add(tracky);
        trackxyz.add(trackz);
        return trackxyz;
    }

    private List<MCParticle> getMCParticlesFromLCIO(EventHeader event){
        List<MCParticle> mcparticles = event.get(MCParticle.class, "MCParticle");
        return mcparticles;
    }

    public int getNSimCalHits(EventHeader event, String simcalhitsCollectionName, MCParticle mcp){

        List<SimCalorimeterHit> simcalhits = event.get(SimCalorimeterHit.class, simcalhitsCollectionName);
        int nHits = 0;
        for(SimCalorimeterHit simcalhit : simcalhits){
            for(int i = 0; i < simcalhit.getMCParticleCount(); i++){
                MCParticle particle = simcalhit.getMCParticle(i);           
                if(particle == mcp)
                    nHits = nHits + 1;
            }
        }
        return nHits;
    }

    private void getClusterMcpMap(List<Cluster> clusters, EventHeader event, boolean verbose, Map<MCParticle, Cluster> mcpClustersMapIn, Map<Cluster, MCParticle> clustersMCPMapIn){

        HashMap<MCParticle, HashSet<Cluster>> mcpClusterMap = new HashMap<MCParticle, HashSet<Cluster>>();
        HashMap<Cluster,MCParticle> mcpClusterMapFinal = new HashMap<Cluster, MCParticle>();

        for(Cluster cluster : clusters){

            double clusterx = cluster.getPosition()[0];
            double clustery = cluster.getPosition()[1];
            double clusterz = cluster.getPosition()[2];
            double clusterEnergy = cluster.getEnergy();

            plots2D.get(String.format("ecal_cluster_positions_xy_plane")).fill(clusterx,clustery);
            plots1D.get(String.format("cluster_truth_stage_0_energy",trackCollectionName)).fill(clusterEnergy);

            //Begin truth matching
            List<MCParticle> mcparticles = getMCParticlesFromLCIO(event);

            List <RawTrackerHit> ecalReadoutHits = event.get(RawTrackerHit.class,ecalReadoutHitsCollectionName);
            Map <MCParticle, int[]>mcParticleMultiplicity = new HashMap<MCParticle, int[]>();
            RelationalTable rawtomc = new BaseRelationalTable(RelationalTable.Mode.MANY_TO_MANY, RelationalTable.Weighting.UNWEIGHTED);
            if (event.hasCollection(LCRelation.class, ecalTruthRelationsName )) {
                List<LCRelation> trueHitRelations = event.get(LCRelation.class, ecalTruthRelationsName);
                for (LCRelation relation : trueHitRelations)
                    if (relation != null && relation.getFrom() != null && relation.getTo() != null){
                        rawtomc.add(relation.getFrom(), relation.getTo());
                    }
                    else System.out.println("failed to build cluster mc relation");
            }

            //Match Cluster->CalorimeterHit(s) to corresponding EcalReadoutHits via
            //CellID match
            List<CalorimeterHit> calorimeterhits = cluster.getCalorimeterHits();
            CalorimeterHit seedhit = ClusterUtilities.findSeedHit(cluster);

            RawTrackerHit readoutMatchHit = null;
            long cellID = seedhit.getCellID();
            for(RawTrackerHit ecalReadoutHit : ecalReadoutHits) {
                long recellID = ecalReadoutHit.getCellID();
                if(cellID == recellID) {
                    readoutMatchHit = ecalReadoutHit;
                }
            }
            if(verbose){
                System.out.println("Checking Cluster with Energy: " + cluster.getEnergy());
                System.out.println("Cluster X: " + cluster.getPosition()[0]);
                System.out.println("Cluster Y: " + cluster.getPosition()[1]);
                System.out.println("Cluster seedhit cellID: " + cellID);
                System.out.println("Cluster seedhit Energy: " + seedhit.getRawEnergy());
            }

            if(readoutMatchHit == null)
                continue;
            if(verbose){
                System.out.println("Matching Readout Hit Found: " + readoutMatchHit.getCellID());
                System.out.println("Readout Hit position: x= " + readoutMatchHit.getPosition()[0] + "; y= " + readoutMatchHit.getPosition()[1] + "; z= " + readoutMatchHit.getPosition()[2]);
            }

            Set<SimCalorimeterHit> simcalhits = rawtomc.allFrom(readoutMatchHit);
            double simcalhit_largestEnergy = 0.0;
            SimCalorimeterHit largest_simcalhit = null;
            double maxMCPEnergy = 0.0;
            for(SimCalorimeterHit simcalhit : simcalhits) {
                if(verbose){
                    System.out.println("Simcalhit energy: " + simcalhit.getRawEnergy());
                }
                if(simcalhit.getRawEnergy() > simcalhit_largestEnergy){
                    simcalhit_largestEnergy = simcalhit.getRawEnergy();
                    largest_simcalhit = simcalhit;
                }
            }
            if(verbose && largest_simcalhit != null){
                System.out.println("Simcalhit with energy: " + largest_simcalhit.getRawEnergy() + " selected");
                System.out.println("Simcalhit cellID: " + largest_simcalhit.getCellID());
                System.out.println("Simcalhit position: x= " + largest_simcalhit.getPosition()[0] + "; y= " + largest_simcalhit.getPosition()[1] + "; z= " + largest_simcalhit.getPosition()[2]);
            }

            double bestMCPEnergyContr = 0.0;
            MCParticle bestMCP = null;

            if(largest_simcalhit == null)
                continue;

            for(int i=0; i < largest_simcalhit.getMCParticleCount(); i++){
                MCParticle mcp = largest_simcalhit.getMCParticle(i);
                double originZ = largest_simcalhit.getMCParticle(i).getOriginZ();
                int PDGID = largest_simcalhit.getMCParticle(i).getPDGID();
                double MCPEnergyFrac = largest_simcalhit.getContributedEnergy(i)/cluster.getEnergy();
                if(verbose){
                    System.out.println("Looping over MCParticles for Simcalhit");
                    System.out.println("MCP energy: " + mcp.getEnergy());
                    System.out.println("MCP energy contribution to simcalhit: " + largest_simcalhit.getContributedEnergy(i));
                    System.out.println("mcp PDGID from mcp = " + mcp.getPDGID());
                    // doesnt work System.out.println("mcp PDGID from simcalhit.getPDG(i) = " + largest_simcalhit.getPDG(i));
                    System.out.println("mcp OriginZ: " + mcp.getOriginZ());
                    System.out.println("mcp EndpointX: " + mcp.getEndPoint().x());
                    System.out.println("mcp EndpointY: " + mcp.getEndPoint().y());
                    System.out.println("mcp EndpointZ: " + mcp.getEndPoint().z());
                    if(mcparticles.contains(mcp)){
                        System.out.println("mcp from simcalhit found in LCIO MCParticle collection");
                    }
                    else
                        System.out.println("mcp from simcalhit NOT FOUND in LCIO MCPartice collection");
                }

                if(mcp.getEnergy() > bestMCPEnergyContr){
                    bestMCPEnergyContr = mcp.getEnergy();
                    bestMCP = mcp;
                }
            }

            if(bestMCP == null)
                continue;

            double distance = Math.sqrt( Math.pow(clusterx - bestMCP.getEndPoint().x(),2) + Math.pow(clustery - bestMCP.getEndPoint().y(),2));
            double energyRatio = cluster.getEnergy()/bestMCP.getEnergy();

            //Cluster has been truth matched to some MCParticle
            plots1D.get(String.format("cluster_truth_stage_1_energy",trackCollectionName)).fill(clusterEnergy);
            plots2D.get("cluster_truth_stage_1_mcpEndpointz_v_ds").fill(bestMCP.getEndPoint().z(),distance);
            plots1D.get(String.format("cluster_truth_stage_1_energy_ratio",this.trackCollectionName)).fill(energyRatio);
            plots2D.get("cluster_truth_stage_1_cluster_v_mcp_energy").fill(cluster.getEnergy(),bestMCP.getEnergy());

            //restrict truth matching to only MCPs that originate at target
            if(bestMCP.getOriginZ() > 0)
                continue;

            plots1D.get(String.format("cluster_truth_stage_2_energy",trackCollectionName)).fill(clusterEnergy);
            plots2D.get("cluster_truth_stage_2_mcpEndpointz_v_ds").fill(bestMCP.getEndPoint().z(),distance);
            plots1D.get(String.format("cluster_truth_stage_2_energy_ratio",this.trackCollectionName)).fill(energyRatio);
            plots2D.get("cluster_truth_stage_2_cluster_v_mcp_energy").fill(cluster.getEnergy(),bestMCP.getEnergy());

            plots2D.get("cluster_truth_stage_2_mcpPy_v_clustery").fill(bestMCP.getPY(),clustery);

            //Require Py sign of MCP to correlate to Top/Bottom cluster
            if(bestMCP.getPY() < 0 & clustery > 0)
                continue;
            if(bestMCP.getPY() > 0 & clustery < 0)
                continue;

            plots1D.get(String.format("cluster_truth_stage_3_energy",trackCollectionName)).fill(clusterEnergy);
            plots2D.get("cluster_truth_stage_3_mcpEndpointz_v_ds").fill(bestMCP.getEndPoint().z(),distance);
            plots1D.get(String.format("cluster_truth_stage_3_energy_ratio",this.trackCollectionName)).fill(energyRatio);
            plots2D.get("cluster_truth_stage_3_cluster_v_mcp_energy").fill(cluster.getEnergy(),bestMCP.getEnergy());

            HashSet<Cluster> c = new HashSet<Cluster>();
            if(mcpClusterMap.containsKey(bestMCP)){
                c = mcpClusterMap.get(bestMCP);
                c.add(cluster);
                mcpClusterMap.put(bestMCP, c);
            }
            else{
                c.add(cluster);
                mcpClusterMap.put(bestMCP, c);
            }
        }

        //Loop over mcpCluster map to check for duplicates
        for(Map.Entry<MCParticle, HashSet<Cluster>> entry : mcpClusterMap.entrySet()){
            MCParticle mcp = entry.getKey();
            HashSet<Cluster> tClusters = entry.getValue();
            plots1D.get("clusters_matched_to_n_mcps").fill(tClusters.size());

            if(tClusters.size() < 2){
                for(Cluster tCluster : tClusters){
                    mcpClusterMapFinal.put(tCluster, mcp);
                    mcpClustersMapIn.put(mcp, tCluster);
                    clustersMCPMapIn.put(tCluster, mcp);
                }
            }

            if(tClusters.size() > 1){
                System.out.println("duplicate clusters matched to single MCP");
                List<Cluster> dupClusters = new ArrayList<Cluster>(tClusters);
                for(int i =0; i < dupClusters.size()-1; i++){
                    plots1D.get("cluster_truth_stage_3_duplicate_mcp_match_dx").fill(dupClusters.get(i).getPosition()[0] - dupClusters.get(i+1).getPosition()[0]);
                    plots1D.get("cluster_truth_stage_3_duplicate_mcp_match_dy").fill(dupClusters.get(i).getPosition()[1] - dupClusters.get(i+1).getPosition()[1]);
                }

                List<Cluster> recoveredDup = new ArrayList<Cluster>();
                for(Cluster tCluster : tClusters){
                    //Cut all matches where MCP.getEndpointZ() is < Ecal face
                    double ecalPosZ = 1300.0;
                    if(mcp.getEndPoint().z() < ecalPosZ)
                        continue;
                    if(tCluster.getPosition()[1] > 0 & mcp.getEndPoint().y() < 0)
                        continue;
                    if(tCluster.getPosition()[1] < 0 & mcp.getEndPoint().y() > 0)
                        continue;
                    if(tCluster.getPosition()[0] > 0 & mcp.getEndPoint().x() < 0)
                        continue;
                    if(tCluster.getPosition()[0] < 0 & mcp.getEndPoint().x() > 0)
                        continue;
                    recoveredDup.add(tCluster);
                    break;
                }

                if(recoveredDup.size() > 1){
                    for(int i =0; i < recoveredDup.size()-1; i++){
                        plots1D.get("cluster_truth_stage_3_cut_remaining_duplicate_mcp_match_dx").fill(recoveredDup.get(i).getPosition()[0] - recoveredDup.get(i+1).getPosition()[0]);
                        plots1D.get("cluster_truth_stage_3_cut_remaining_duplicate_mcp_match_dy").fill(recoveredDup.get(i).getPosition()[1] - recoveredDup.get(i+1).getPosition()[1]);

                    }
                }
                else if(recoveredDup.size() > 0 & recoveredDup.size() < 2) {
                    mcpClusterMapFinal.put(recoveredDup.get(0),mcp);
                    mcpClustersMapIn.put(mcp, recoveredDup.get(0));
                    clustersMCPMapIn.put(recoveredDup.get(0),mcp);
                }
            }
        }

        for(Map.Entry<Cluster, MCParticle> entry: mcpClusterMapFinal.entrySet()){
            Cluster cluster = entry.getKey();
            MCParticle mcp = entry.getValue();
            double distance = Math.sqrt( Math.pow(cluster.getPosition()[0] - mcp.getEndPoint().x(),2) + Math.pow(cluster.getPosition()[1] - mcp.getEndPoint().y(),2));
            double energyRatio = cluster.getEnergy()/mcp.getEnergy();
            plots1D.get(String.format("cluster_truth_stage_final_energy",trackCollectionName)).fill(cluster.getEnergy());
            plots2D.get("cluster_truth_stage_final_mcpEndpointz_v_ds").fill(mcp.getEndPoint().z(),distance);
            plots1D.get(String.format("cluster_truth_stage_final_energy_ratio",this.trackCollectionName)).fill(energyRatio);
            plots2D.get("cluster_truth_stage_final_cluster_v_mcp_energy").fill(cluster.getEnergy(),mcp.getEnergy());
        }
        //return mcpClusterMapFinal;
    }

    public int getNSimTrackerHits(EventHeader event, MCParticle mcp){
        //Check how many hits this MCP left in the tracker
        int nmcpHits = 0;
        Set<Integer> layers = new HashSet<Integer>();
        List<SimTrackerHit> simhits = event.get(SimTrackerHit.class, "TrackerHits");
        for(SimTrackerHit simhit : simhits){
            if(layers.contains(simhit.getLayerNumber()))
                continue;
            MCParticle simhitmcp = simhit.getMCParticle();
            if(simhitmcp == mcp){
                layers.add(simhit.getLayerNumber());
                nmcpHits = nmcpHits + 1;
            }
        }
        return nmcpHits;
    }

    public Track getMcpBestTrack(MCParticle mcp, Map<MCParticle, Map<Track, Integer>> mcpTrackMap){
        int mostHits = 0;
        Track mcpBestTrack = null;
        for(Map.Entry<Track, Integer> entry : mcpTrackMap.get(mcp).entrySet()){
            if(entry.getValue() > mostHits){
                mostHits = entry.getValue();
                mcpBestTrack = entry.getKey();
            }
        }
        return mcpBestTrack;
    }
}


