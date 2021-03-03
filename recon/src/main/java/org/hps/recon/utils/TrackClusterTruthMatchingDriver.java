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

import org.hps.recon.tracking.TrackUtils;
import org.hps.recon.tracking.TrackData;
import org.hps.record.StandardCuts;
import org.hps.recon.ecal.cluster.ClusterUtilities;
//import org.hps.recon.ecal.cluster.ClusterCorrectionUtilities;

import org.hps.util.Pair;
import org.hps.util.RK4integrator;

import org.lcsim.event.EventHeader;
import org.lcsim.event.LCRelation;
import org.lcsim.event.MCParticle;
import org.lcsim.event.SimTrackerHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.Track;
import org.lcsim.event.TrackerHit;
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

    private org.lcsim.geometry.FieldMap fM;
    private ITree tree;
    private IHistogramFactory histogramFactory;
    private Map<String, IHistogram1D> plots1D;
    private Map<String, IHistogram2D> plots2D;
    String[] identifiers = {"truth_matched"};

    RelationalTable hitToRotated = null;
    RelationalTable hitToStrips = null;
    RelationalTable TrktoData = new BaseRelationalTable(RelationalTable.Mode.ONE_TO_ONE, RelationalTable.Weighting.UNWEIGHTED);
    HPSEcal3 ecal;
    double beamEnergy = 2.3;
    protected StandardCuts cuts = new StandardCuts();
    TrackClusterMatcherInter matcher;

    //INITIALIZE TERMS TO EVALUATE MATCHING PERFORMANCE
    int ntruthpairs_ele = 0;
    int ntruthpairs_pos = 0;
    int ntruthTracks_ele = 0;
    int ntruthTracks_pos = 0;
    int nrecoTracks_ele = 0;
    int nrecoTracks_pos = 0;
    int ntruthClusters = 0;
    int nClusters = 0;

    double goodMatch_ele=0;
    double fakeMatch_ele=0;
    double noMatch_ele=0;
    double missedMatch_ele=0;
    double unknownMatch_ele=0;
    double goodMatch_pos=0;
    double fakeMatch_pos=0;
    double noMatch_pos=0;
    double missedMatch_pos=0;
    double unknownMatch_pos=0;

    double noMatchTruthTrack_ele=0;
    double unknownMatchTruthTrack_ele=0;
    double unknownMatchTruthCluster_ele=0;
    double truthTrackOutsideEcal_ele=0; 
    double trackOutsideEcal_ele=0;
    double noMatchTruthTrack_pos=0;
    double unknownMatchTruthTrack_pos=0;
    double unknownMatchTruthCluster_pos=0;
    double truthTrackOutsideEcal_pos=0; 
    double trackOutsideEcal_pos=0;

    double efficiency_ele;
    double efficiency_pos;
    double fakerate_ele;
    double fakerate_pos;


    //Counting fake rate
    boolean verbose = true;
    double NtruthEleClustPairs = 0;
    double NtruthPosClustPairs = 0;

    //Check track efficiency 
    double NpossibleTracks = 0;
    double NrecoTruthTracks = 0;
    double trackEfficiency = 0;

    //check number of reco tracks that are matched to the same MCP "duplicates"
    double nDuplicates = 0;

    //time difference for track and cluster
    double trackClusterTimeOffset;
    boolean truthComparisons = true;


    //Collection Names
    String ecalScoringPlaneHitsCollectionName = "TrackerHitsECal";
    String trackCollectionName = "KalmanFullTracks";
    String trackToScoringPlaneHitRelationsName = "TrackToEcalScoringPlaneHitRelations";

    String ecalClustersCollectionName = "EcalClusters";
    String ecalReadoutHitsCollectionName = "EcalReadoutHits";
    String ecalTruthRelationsName = "EcalTruthRelations";

    private Set<SimTrackerHit> simhitsontrack = new HashSet<SimTrackerHit>();


    public void setTrackClusterTimeOffset(double input) {
        trackClusterTimeOffset = input;
    }
    public void setTrackCollectionName(String trackCollectionName) {
        this.trackCollectionName = trackCollectionName;
    }

    public void saveHistograms() {
        System.out.println("[TrackClusterTruthMatchingDriver] Saving Histogram for " + this.trackCollectionName);
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




//FINAL CLUSTER TRUTH MATCHING PLOTS
        
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


        //duplicate matches
        plots1D.put(String.format("cluster_truth_stage_3_duplicate_mcp_match_dx",trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_truth_stage_3_duplicate_mcp_match_dx",trackCollectionName),  1000, -1000, 1000));
        plots1D.put(String.format("cluster_truth_stage_3_duplicate_mcp_match_dy",trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_truth_stage_3_duplicate_mcp_match_dy",trackCollectionName),  1000, -1000, 1000));

        plots1D.put(String.format("cluster_truth_stage_3_cut_remaining_duplicate_mcp_match_dx",trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_truth_stage_3_cut_remaining_duplicate_mcp_match_dx",trackCollectionName),  1000, -1000, 1000));
        plots1D.put(String.format("cluster_truth_stage_3_cut_remaining_duplicate_mcp_match_dy",trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_truth_stage_3_cut_remaining_duplicate_mcp_match_dy",trackCollectionName),  1000, -1000, 1000));

        plots1D.put(String.format("cluster_truth_stage_final_energy",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_truth_stage_final_energy",this.trackCollectionName), 1000, 0, 10));
        plots2D.put(String.format( "cluster_truth_stage_final_mcpEndpointz_v_ds"), histogramFactory.createHistogram2D(String.format("cluster_truth_stage_final_mcpEndpointz_v_ds"),200, -100, 1900,1000, 0, 2000));
        plots1D.put(String.format("cluster_truth_stage_final_energy_ratio",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_truth_stage_final_energy_ratio",this.trackCollectionName), 1000, 0, 10));
        plots2D.put(String.format("cluster_truth_stage_final_cluster_v_mcp_energy",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("cluster_truth_stage_final_cluster_v_mcp_energy",this.trackCollectionName),1000, 0, 5,1000, 0, 5));

        plots1D.put(String.format("nclusters_matched_to_single_mcp"), histogramFactory.createHistogram1D(String.format("nclusters_matched_to_single_mcp"),  10, 0, 10));



/////////////////////////////

//PLOTS TO CHECK MCPARTICLES THAT ARE MATCHED TO MULTIPLE CLUSTERS...THERE ARE
//CASES WHERE TWO CLUSTERS IN OPPOSITE SIDES OF THE ECAL ARE MATCHED TO THE
//SAME MCP....SHOULD BE IMPOSSIBLE


        plots2D.put(String.format("cluster_mcp_isBackscatter_zy"), histogramFactory.createHistogram2D(String.format("cluster_mcp_isBackscatter_zy"),1000, -100, 1900,1000, -500, 500));


        plots2D.put(String.format("clustery_v_mcpy",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("clustery_v_mcpy",this.trackCollectionName),1000, -500, 500,1000, -500, 500));
        plots2D.put(String.format("clusterx_v_mcpx",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("clusterx_v_mcpx",this.trackCollectionName),1000, -500, 500,1000, -500, 500));

        plots2D.put(String.format("cluster_vs_photon_mcp_truth_energy",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("cluster_vs_photon_mcp_truth_energy",this.trackCollectionName),1000, 0, 5,1000, 0, 5));
        plots2D.put(String.format("cluster_vs_electron_mcp_truth_energy",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("cluster_vs_electron_mcp_truth_energy",this.trackCollectionName),1000, 0, 5,1000, 0, 5));
        plots2D.put(String.format("cluster_vs_positron_mcp_truth_energy",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("cluster_vs_positron_mcp_truth_energy",this.trackCollectionName),1000, 0, 5,1000, 0, 5));


        plots2D.put(String.format("cluster_vs_photon_mcp_truth_energy_onEdge",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("cluster_vs_photon_mcp_truth_energy_onEdge",this.trackCollectionName),1000, 0, 5,1000, 0, 5));
        plots2D.put(String.format("cluster_vs_electron_mcp_truth_energy_onEdge",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("cluster_vs_electron_mcp_truth_energy_onEdge",this.trackCollectionName),1000, 0, 5,1000, 0, 5));
        plots2D.put(String.format("cluster_vs_positron_mcp_truth_energy_onEdge",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("cluster_vs_positron_mcp_truth_energy_onEdge",this.trackCollectionName),1000, 0, 5,1000, 0, 5));

        plots2D.put(String.format("cluster_vs_photon_mcp_truth_energy_offEdge",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("cluster_vs_photon_mcp_truth_energy_offEdge",this.trackCollectionName),1000, 0, 5,1000, 0, 5));
        plots2D.put(String.format("cluster_vs_electron_mcp_truth_energy_offEdge",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("cluster_vs_electron_mcp_truth_energy_offEdge",this.trackCollectionName),1000, 0, 5,1000, 0, 5));
        plots2D.put(String.format("cluster_vs_positron_mcp_truth_energy_offEdge",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("cluster_vs_positron_mcp_truth_energy_offEdge",this.trackCollectionName),1000, 0, 5,1000, 0, 5));

        plots2D.put(String.format("cluster_vs_photon_mcp_truth_energy_unsureEdge",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("cluster_vs_photon_mcp_truth_energy_unsureEdge",this.trackCollectionName),1000, 0, 5,1000, 0, 5));
        plots2D.put(String.format("cluster_vs_electron_mcp_truth_energy_unsureEdge",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("cluster_vs_electron_mcp_truth_energy_unsureEdge",this.trackCollectionName),1000, 0, 5,1000, 0, 5));
        plots2D.put(String.format("cluster_vs_positron_mcp_truth_energy_unsureEdge",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("cluster_vs_positron_mcp_truth_energy_unsureEdge",this.trackCollectionName),1000, 0, 5,1000, 0, 5));


        plots2D.put(String.format("cluster_positron_mcp_duplicates_dx_v_dy"), histogramFactory.createHistogram2D(String.format("cluster_positron_mcp_duplicates_dx_v_dy"),700, -700, 700, 700, -700, 700));
        plots2D.put(String.format("cluster_photon_mcp_duplicates_dx_v_dy"), histogramFactory.createHistogram2D(String.format("cluster_photon_mcp_duplicates_dx_v_dy"),700, -700, 700, 700, -700, 700));
        plots2D.put(String.format("cluster_electron_mcp_duplicates_dx_v_dy"), histogramFactory.createHistogram2D(String.format("cluster_electron_mcp_duplicates_dx_v_dy"),700, -700, 700, 700, -700, 700));

        plots1D.put(String.format("cluster_photon_mcp_truth_originz"), histogramFactory.createHistogram1D(String.format("cluster_photon_mcp_truth_originz"),  1100, -200, 2000));
        plots1D.put(String.format("cluster_electron_mcp_truth_originz"), histogramFactory.createHistogram1D(String.format("cluster_electron_mcp_truth_originz"),  1100, -200, 2000));
        plots1D.put(String.format("cluster_positron_mcp_truth_originz"), histogramFactory.createHistogram1D(String.format("cluster_positron_mcp_truth_originz"),  1100, -200, 2000));

        plots1D.put(String.format("cluster_photon_mcp_truth_energy"), histogramFactory.createHistogram1D(String.format("cluster_photon_mcp_truth_energy"),  500, 0, 5));
        plots1D.put(String.format("cluster_electron_mcp_truth_energy"), histogramFactory.createHistogram1D(String.format("cluster_electron_mcp_truth_energy"),  500, 0, 5));
        plots1D.put(String.format("cluster_positron_mcp_truth_energy"), histogramFactory.createHistogram1D(String.format("cluster_positron_mcp_truth_energy"),  500, 0, 5));

        plots1D.put(String.format("cluster_mcp_truth_match_event_multiplicity"), histogramFactory.createHistogram1D(String.format("cluster_mcp_truth_match_event_multiplicity"),  20, 0, 20));
        plots1D.put(String.format("cluster_event_multiplicity"), histogramFactory.createHistogram1D(String.format("cluster_event_multiplicity"),  20, 0, 20));



        plots1D.put(String.format("cluster_photon_mcp_duplicates_dx",trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_photon_mcp_duplicates_dx",trackCollectionName),  1000, -1000, 1000));
        plots1D.put(String.format("cluster_photon_mcp_duplicates_dy",trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_photon_mcp_duplicates_dy",trackCollectionName),  1000, -1000, 1000));
        plots1D.put(String.format("cluster_photon_mcp_duplicates_dr",trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_photon_mcp_duplicates_dr",trackCollectionName),  1000, -1000, 1000));
        plots1D.put(String.format("cluster_photon_mcp_duplicates_dE",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_photon_mcp_duplicates_dE",this.trackCollectionName), 1000, 0, 10));
        plots1D.put(String.format("cluster_photon_mcp_duplicates_originz"), histogramFactory.createHistogram1D(String.format("cluster_photon_mcp_duplicates_originz"),  1100, -200, 2000));
        plots2D.put(String.format("cluster_photon_mcp_duplicates_xy"), histogramFactory.createHistogram2D(String.format("cluster_photon_mcp_duplicates_xy"),1000, -500, 500,1000, -500, 500));

        plots1D.put(String.format("cluster_electron_mcp_duplicates_dx",trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_electron_mcp_duplicates_dx",trackCollectionName),  1000, -1000, 1000));
        plots1D.put(String.format("cluster_electron_mcp_duplicates_dy",trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_electron_mcp_duplicates_dy",trackCollectionName),  1000, -1000, 1000));
        plots1D.put(String.format("cluster_electron_mcp_duplicates_dr",trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_electron_mcp_duplicates_dr",trackCollectionName),  1000, -1000, 1000));
        plots1D.put(String.format("cluster_electron_mcp_duplicates_dE",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_electron_mcp_duplicates_dE",this.trackCollectionName), 1000, 0, 10));
        plots1D.put(String.format("cluster_electron_mcp_duplicates_originz"), histogramFactory.createHistogram1D(String.format("cluster_electron_mcp_duplicates_originz"),  1100, -200, 2000));
        plots2D.put(String.format("cluster_electron_mcp_duplicates_xy"), histogramFactory.createHistogram2D(String.format("cluster_electron_mcp_duplicates_xy"),1000, -500, 500,1000, -500, 500));

        plots1D.put(String.format("cluster_positron_mcp_duplicates_dx",trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_positron_mcp_duplicates_dx",trackCollectionName),  1000, -1000, 1000));
        plots1D.put(String.format("cluster_positron_mcp_duplicates_dy",trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_positron_mcp_duplicates_dy",trackCollectionName),  1000, -1000, 1000));
        plots1D.put(String.format("cluster_positron_mcp_duplicates_dr",trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_positron_mcp_duplicates_dr",trackCollectionName),  1000, -1000, 1000));
        plots1D.put(String.format("cluster_positron_mcp_duplicates_dE",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_positron_mcp_duplicates_dE",this.trackCollectionName), 1000, 0, 10));
        plots1D.put(String.format("cluster_positron_mcp_duplicates_originz"), histogramFactory.createHistogram1D(String.format("cluster_positron_mcp_duplicates_originz"),  1100, -200, 2000));
        plots2D.put(String.format("cluster_positron_mcp_duplicates_xy"), histogramFactory.createHistogram2D(String.format("cluster_positron_mcp_duplicates_xy"),1000, -500, 500,1000, -500, 500));

///////////////////////////////////////////////////////////////////////////////////


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

        //Track reconstruction efficiency
        plots1D.put(String.format("%s_track_reco_efficiency",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_track_reco_efficiency",this.trackCollectionName), 10, 0, 10));
        //Number of duplicate MCPs matched to reco Tracks
        plots1D.put(String.format("%s_n_duplicate_TrackMCP_matches",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_n_duplicate_TrackMCP_matches",this.trackCollectionName), 10, 0, 10));

        //Checking the track-cluster matching fake rate and matching efficiency
        plots1D.put(String.format("%s_ele_fakeRate",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_fakeRate",this.trackCollectionName), 10, 0, 10));
        plots1D.put(String.format("%s_pos_fakeRate",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_fakeRate",this.trackCollectionName), 10, 0, 10));
        plots1D.put(String.format("%s_ele_Efficiency",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_Efficiency",this.trackCollectionName), 10, 0, 10));
        plots1D.put(String.format("%s_pos_Efficiency",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_Efficiency",this.trackCollectionName), 10, 0, 10));


//TRUTH MATCHED TRACKS WITH MCPARTICLES

        plots2D.put(String.format("%s_ele_mcp_v_truth_track_momentum",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_mcp_v_truth_track_momentum",this.trackCollectionName),600, -1, 5, 600, -1, 5));
        plots2D.put(String.format("%s_pos_mcp_v_truth_track_momentum",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_mcp_v_truth_track_momentum",this.trackCollectionName),600, -1, 5, 600, -1, 5));

//PLOTS FOR TRUTH TRACK CLUSTER PAIRS WITH DIFFERENT ACCEPTANCES

        //TRUTH TRACK CLUSTER PAIRS
        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_momentum",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_momentum",this.trackCollectionName), 1000, 0, 5));
            //momentum
        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_momentum",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_momentum",this.trackCollectionName), 1000, 0, 5));

        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_px",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_px",this.trackCollectionName), 600, -1, 5));
        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_py",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_py",this.trackCollectionName), 600, -1, 5));
        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_pz",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_pz",this.trackCollectionName), 600, -1, 5));

        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_px",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_px",this.trackCollectionName), 600, -1, 5));
        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_py",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_py",this.trackCollectionName), 600, -1, 5));
        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_pz",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_pz",this.trackCollectionName), 600, -1, 5));

        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_dx",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_dx",trackCollectionName),  800, -200, 200));
        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_dy",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_dy",trackCollectionName),  800, -200, 200));
        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_dx",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_dx",trackCollectionName),  800, -200, 200));
        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_dy",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_dy",trackCollectionName),  800, -200, 200));
        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_dz",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_dz",trackCollectionName),  800, -200, 200));
        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_dz",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_dz",trackCollectionName),  800, -200, 200));
        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_EdivP",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_EdivP",trackCollectionName),  1000, 0, 10));
        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_EdivP",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_EdivP",trackCollectionName),  1000, 0, 10));
        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_clusterMCP_Eratio",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_clusterMCP_Eratio",this.trackCollectionName), 1000, 0, 10));
        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_clusterMCP_Eratio",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_clusterMCP_Eratio",this.trackCollectionName), 1000, 0, 10));
        plots2D.put(String.format("%s_ele_truth_track_cluster_pair_xypos",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_truth_track_cluster_pair_xypos",this.trackCollectionName),1000, -500, 500,1000, -500, 500));
        plots2D.put(String.format("%s_pos_truth_track_cluster_pair_xypos",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_truth_track_cluster_pair_xypos",this.trackCollectionName),1000, -500, 500,1000, -500, 500));

        //TEST FOR CHECKING TRACKS MATCHED TO MULTIPLE TRUTH CLUSTERS
        //testing_track_w_multiple_cluster_matches_dx
        plots1D.put(String.format("testing_track_w_multiple_cluster_matches_dx",trackCollectionName), histogramFactory.createHistogram1D(String.format("testing_track_w_multiple_cluster_matches_dx",trackCollectionName),  800, -200, 200));
        plots1D.put(String.format("testing_track_w_multiple_cluster_matches_dy",trackCollectionName), histogramFactory.createHistogram1D(String.format("testing_track_w_multiple_cluster_matches_dy",trackCollectionName),  800, -200, 200));
        plots1D.put(String.format("testing_track_w_multiple_cluster_matches_dz",trackCollectionName), histogramFactory.createHistogram1D(String.format("testing_track_w_multiple_cluster_matches_dz",trackCollectionName),  800, -200, 200));
        plots1D.put(String.format("testing_track_w_multiple_cluster_matches_dt",trackCollectionName), histogramFactory.createHistogram1D(String.format("testing_track_w_multiple_cluster_matches_dt",trackCollectionName),  800, -200, 200));
        plots1D.put(String.format("testing_track_w_multiple_cluster_matches_dr",trackCollectionName), histogramFactory.createHistogram1D(String.format("testing_track_w_multiple_cluster_matches_dr",trackCollectionName),  800, -200, 200));

        plots1D.put(String.format("testing_track_w_multiple_cluster_matches_EdivP",trackCollectionName), histogramFactory.createHistogram1D(String.format("testing_track_w_multiple_cluster_matches_EdivP",trackCollectionName),  1000, 0, 10));

        plots1D.put(String.format("testing_track_w_multiple_cluster_matches_clusterMCP_Eratio",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("testing_track_w_multiple_cluster_matches_clusterMCP_Eratio",this.trackCollectionName), 1000, 0, 10));
        //END

        //TRUTH TRACK CLUSTER PAIRS INSIDE ECAL
        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_insideEcal_momentum",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_insideEcal_momentum",this.trackCollectionName), 1000, 0, 5));
        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_insideEcal_momentum",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_insideEcal_momentum",this.trackCollectionName), 1000, 0, 5));


        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_insideEcal_dx",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_insideEcal_dx",trackCollectionName),  800, -200, 200));
        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_insideEcal_dy",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_insideEcal_dy",trackCollectionName),  800, -200, 200));
        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_insideEcal_dx",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_insideEcal_dx",trackCollectionName),  800, -200, 200));
        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_insideEcal_dy",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_insideEcal_dy",trackCollectionName),  800, -200, 200));
        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_insideEcal_dz",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_insideEcal_dz",trackCollectionName),  800, -200, 200));
        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_insideEcal_dz",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_insideEcal_dz",trackCollectionName),  800, -200, 200));
        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_insideEcal_EdivP",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_insideEcal_EdivP",trackCollectionName),  1000, 0, 10));
        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_insideEcal_EdivP",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_insideEcal_EdivP",trackCollectionName),  1000, 0, 10));
        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_insideEcal_clusterMCP_Eratio",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_insideEcal_clusterMCP_Eratio",this.trackCollectionName), 1000, 0, 10));
        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_insideEcal_clusterMCP_Eratio",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_insideEcal_clusterMCP_Eratio",this.trackCollectionName), 1000, 0, 10));

            //cluster and track xy plots
            
        plots2D.put(String.format("%s_ele_truth_track_cluster_pair_insideEcal_xypos",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_truth_track_cluster_pair_insideEcal_xypos",this.trackCollectionName),1000, -500, 500,1000, -500, 500));
        plots2D.put(String.format("%s_pos_truth_track_cluster_pair_insideEcal_xypos",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_truth_track_cluster_pair_insideEcal_xypos",this.trackCollectionName),1000, -500, 500,1000, -500, 500));
        plots2D.put(String.format("%s_ele_truth_track_cluster_pair_insideEcal_cluster_xypos",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_truth_track_cluster_pair_insideEcal_cluster_xypos",this.trackCollectionName),1000, -500, 500,1000, -500, 500));
        plots2D.put(String.format("%s_pos_truth_track_cluster_pair_insideEcal_cluster_xypos",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_truth_track_cluster_pair_insideEcal_cluster_xypos",this.trackCollectionName),1000, -500, 500,1000, -500, 500));

        //TRUTH TRACK CLUSTER PAIRS OUTSIDE ECAL
        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_momentum",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_momentum",this.trackCollectionName), 1000, 0, 5));
        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_momentum",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_momentum",this.trackCollectionName), 1000, 0, 5));
        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_Pz",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_Pz",this.trackCollectionName), 600, -1, 5));
        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_Pz",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_Pz",this.trackCollectionName), 600, -1, 5));

        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_dx",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_dx",trackCollectionName),  800, -200, 200));
        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_dy",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_dy",trackCollectionName),  800, -200, 200));
        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_dx",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_dx",trackCollectionName),  800, -200, 200));
        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_dy",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_dy",trackCollectionName),  800, -200, 200));
        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_dz",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_dz",trackCollectionName),  800, -200, 200));
        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_dz",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_dz",trackCollectionName),  800, -200, 200));
        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_EdivP",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_EdivP",trackCollectionName),  1000, 0, 10));
        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_EdivP",trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_EdivP",trackCollectionName),  1000, 0, 10));
        plots1D.put(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_clusterMCP_Eratio",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_clusterMCP_Eratio",this.trackCollectionName), 1000, 0, 10));
        plots1D.put(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_clusterMCP_Eratio",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_clusterMCP_Eratio",this.trackCollectionName), 1000, 0, 10));
        plots2D.put(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_xypos",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_xypos",this.trackCollectionName),1000, -500, 500,1000, -500, 500));
        plots2D.put(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_xypos",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_xypos",this.trackCollectionName),1000, -500, 500,1000, -500, 500));
        plots2D.put(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_cluster_xypos",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_cluster_xypos",this.trackCollectionName),1000, -500, 500,1000, -500, 500));
        plots2D.put(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_cluster_xypos",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_cluster_xypos",this.trackCollectionName),1000, -500, 500,1000, -500, 500));


        //track and truthtrack momentum
        plots1D.put(String.format("%s_ele_track_momentum",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_track_momentum",this.trackCollectionName), 1000, 0, 5));
        plots1D.put(String.format("%s_pos_track_momentum",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_track_momentum",this.trackCollectionName), 1000, 0, 5));

        plots1D.put(String.format("%s_ele_delta_momentum",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_delta_momentum",this.trackCollectionName), 400, 0, 4));
        plots1D.put(String.format("%s_pos_delta_momentum",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_delta_momentum",this.trackCollectionName), 400, 0, 4));

        plots1D.put(String.format("%s_ele_truth_track_momentum",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_momentum",this.trackCollectionName), 1000, 0, 5));
        plots1D.put(String.format("%s_pos_truth_track_momentum",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_momentum",this.trackCollectionName), 1000, 0, 5));
        plots1D.put(String.format("%s_pos_truth_track_MCP_momentum_ratio",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_truth_track_MCP_momentum_ratio",this.trackCollectionName), 1000, 0, 5));
        plots1D.put(String.format("%s_ele_truth_track_MCP_momentum_ratio",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_truth_track_MCP_momentum_ratio",this.trackCollectionName), 1000, 0, 5));


        //duplicate track momentum
        plots1D.put(String.format("%s_ele_duplicate_track_momentum",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_duplicate_track_momentum",this.trackCollectionName), 1000, 0, 5));
        plots1D.put(String.format("%s_pos_duplicate_track_momentum",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_duplicate_track_momentum",this.trackCollectionName), 1000, 0, 5));
        //Checking E/P for track+cluster pair
        
        plots1D.put(String.format("cluster_truth_energy",trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_truth_energy",trackCollectionName),  100, 0, 5));

        //MCP hits along z
        plots1D.put(String.format("%s_EcalCluster_Simcalhit_MCParticle_origin_z",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_EcalCluster_Simcalhit_MCParticle_origin_z",this.trackCollectionName), 1500, -1000, 2000));

        //Hit multiplicity for truth matching
        plots1D.put(String.format("%s_ele_track_maxMCPmultiplicity",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_track_maxMCPmultiplicity",this.trackCollectionName),16, 0, 16));
        plots1D.put(String.format("%s_pos_track_maxMCPmultiplicity",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_track_maxMCPmultiplicity",this.trackCollectionName),16, 0, 16));

        //track positions at Ecal in XY plane
        plots2D.put(String.format("%s_ele_track_xypos_at_ecal",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_track_xypos_at_ecal",this.trackCollectionName),1000, -500, 500,1000, -500, 500));
        plots2D.put(String.format("%s_pos_track_xypos_at_ecal",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_track_xypos_at_ecal",this.trackCollectionName),1000, -500, 500,1000, -500, 500));

        plots2D.put(String.format("%s_ele_truth_track_xypos_at_ecal",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_ele_truth_track_xypos_at_ecal",this.trackCollectionName),1000, -500, 500,1000, -500, 500));
        plots2D.put(String.format("%s_pos_truth_track_xypos_at_ecal",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("%s_pos_truth_track_xypos_at_ecal",this.trackCollectionName),1000, -500, 500,1000, -500, 500));



        //Cluster Positions in XY plane
        plots2D.put(String.format("ecal_cluster_positions_xy_plane",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("ecal_cluster_positions_xy_plane",this.trackCollectionName),1000, -500, 500,1000, -500, 500));
        plots2D.put(String.format("cluster_truthpositions_xy_plane",this.trackCollectionName), histogramFactory.createHistogram2D(String.format("cluster_truthpositions_xy_plane",this.trackCollectionName),50, -280, 370, 18, -110, 124));

        //Check track quality
        plots1D.put(String.format("%s_ele_track_chi2divndf",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_ele_track_chi2divndf",this.trackCollectionName), 200, 0, 200));
        plots1D.put(String.format("%s_pos_track_chi2divndf",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("%s_pos_track_chi2divndf",this.trackCollectionName), 200, 0, 200));

        plots1D.put(String.format("cluster_truthMCP_energy_ratio_sorted",this.trackCollectionName), histogramFactory.createHistogram1D(String.format("cluster_truthMCP_energy_ratio_sorted",this.trackCollectionName), 1000, 0, 10));



    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    int localEventNumber = 0;
    public void startOfData() {
        System.out.println("Starting job");
        bookHistograms();

        //init matcher being tested
        matcher = TrackClusterMatcherFactory.create("TrackClusterMatcher2019");
        matcher.setTrackCollectionName(this.trackCollectionName);
        matcher.enablePlots(true);
    }

    public void endOfData() {

        matcher.saveHistograms();
       
        //track reco efficiency
        trackEfficiency = NrecoTruthTracks/NpossibleTracks;
        for(int i=0; i < Math.round(trackEfficiency*100); i++)
            plots1D.get(String.format("%s_track_reco_efficiency",this.trackCollectionName)).fill(1.0);
        System.out.println("LOOK! nDuplicates = " + nDuplicates);
        for(int i = 0; i < Math.round(nDuplicates); i++)
            plots1D.get(String.format("%s_n_duplicate_TrackMCP_matches",this.trackCollectionName)).fill(1.0);

        /*
        for(int i=0; i < Math.round(fakerate_ele*100); i++){
            System.out.println("Filling Histogram fakerate");
            plots1D.get(String.format("%s_ele_fakeRate",this.trackCollectionName)).fill(1.0);
        }
        for(int i=0; i < Math.round(fakerate_pos*100); i++)
            plots1D.get(String.format("%s_pos_fakeRate",this.trackCollectionName)).fill(1.0);
        for(int i=0; i < Math.round(efficiency_ele*100); i++)
            plots1D.get(String.format("%s_pos_Efficiency",this.trackCollectionName)).fill(1.0);
        for(int i=0; i < Math.round(efficiency_pos*100); i++)
            plots1D.get(String.format("%s_ele_Efficiency",this.trackCollectionName)).fill(1.0);
            */
        saveHistograms();
        System.out.println("goodMatch_ele: " + goodMatch_ele);
        System.out.println("fakeMatch_ele: " + fakeMatch_ele);
        System.out.println("noMatch_ele: " + noMatch_ele);
        System.out.println("missedMatch_ele: " + missedMatch_ele);
        System.out.println("unknownMatch_ele: " + unknownMatch_ele);
        System.out.println("goodMatch_pos: " + goodMatch_pos);
        System.out.println("fakeMatch_pos: " + fakeMatch_pos);
        System.out.println("noMatch_pos: " + noMatch_pos);
        System.out.println("missedMatch_pos: " + missedMatch_pos);
        System.out.println("unknownMatch_pos: " + unknownMatch_pos);

        System.out.println("noMatchTruthTrack_ele: " + noMatchTruthTrack_ele);
        System.out.println("unknownMatchTruthTrack_ele: " + unknownMatchTruthTrack_ele);
        System.out.println("unknownMatchTruthCluster_ele: " + unknownMatchTruthCluster_ele);
        System.out.println("truthTrackOutsideEcal_ele: " + truthTrackOutsideEcal_ele);
        System.out.println("trackOutsideEcal_ele: " + trackOutsideEcal_ele);
        System.out.println("noMatchTruthTrack_pos: " + noMatchTruthTrack_pos);
        System.out.println("unknownMatchTruthTrack_pos: " + unknownMatchTruthTrack_pos);
        System.out.println("unknownMatchTruthCluster_pos: " + unknownMatchTruthCluster_pos);
        System.out.println("truthTrackOutsideEcal_pos: " + truthTrackOutsideEcal_pos);
        System.out.println("trackOutsideEcal_pos: " + trackOutsideEcal_pos);
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

    public List<MCParticle> getPossibleTrackMCPs(EventHeader event, int minHitsOnTrack){
        //Loop over all simhits in the LCIO file and check how many Tracks
        //could possibly be reconstructed...where a possible track is counted
        //a MC particle leaves minHitsOnTrack hits
        List<SimTrackerHit> simhits =  event.get(SimTrackerHit.class, "TrackerHits");
        Map<MCParticle, int[]> nMCPhits = new HashMap<MCParticle, int[]>();
        for(SimTrackerHit simhit : simhits){
            MCParticle mcp = simhit.getMCParticle();
            if(!nMCPhits.containsKey(mcp)){
                nMCPhits.put(mcp, new int[1]);
                nMCPhits.get(mcp)[0] = 0;
            }
            nMCPhits.get(mcp)[0]++;
        }

        List<MCParticle> possibleTracks = new ArrayList<MCParticle>();
        for(Map.Entry<MCParticle, int[]> entry : nMCPhits.entrySet()){
            if(entry.getValue()[0] > minHitsOnTrack ){
                possibleTracks.add(entry.getKey());
            }
        }

        return possibleTracks;
    }

    protected void process(EventHeader event) {

        this.localEventNumber = this.localEventNumber + 1;
        cuts.setTrackClusterTimeOffset(44.8);

        //If event has no collection of tracks, skip
        if(!event.hasCollection(Track.class, trackCollectionName)) return;

        //If even doesnt have collection of Ecal scoring plane hits, skip
        if(!event.hasCollection(SimTrackerHit.class, ecalScoringPlaneHitsCollectionName)) return;

        //Get EcalClusters from event
        List<Cluster> clusters = event.get(Cluster.class, ecalClustersCollectionName);

        //Truth Match Clusters to MCParticles    
        //Map<Cluster,MCParticle> truthClustersMap = getUniqueTruthClusters(event,clusters);
        Map<Cluster, MCParticle> truthClustersMap = getClusterMcpMap(clusters, event, false);
        List<Cluster> truthClusters = new ArrayList<Cluster>();
        for(Map.Entry<Cluster,MCParticle> entry : truthClustersMap.entrySet()){
            truthClusters.add(entry.getKey());
        }
        drawEcalFace(truthClusters);

        for(Cluster cluster : truthClusters){
            double clustx = cluster.getPosition()[0];
            double clusty = cluster.getPosition()[1];
            double clustz = cluster.getPosition()[2];
            double clusterEnergy = cluster.getEnergy();
            double energyRatio = clusterEnergy/truthClustersMap.get(cluster).getEnergy();

            plots1D.get(String.format("cluster_truth_energy",trackCollectionName)).fill(clusterEnergy);
            plots2D.get(String.format("cluster_truthpositions_xy_plane")).fill(clustx,clusty);
            plots1D.get(String.format("cluster_truthMCP_energy_ratio_sorted",this.trackCollectionName)).fill(energyRatio);
        }


        // Get collection of tracks from event
        List<Track> tracks = event.get(Track.class, trackCollectionName);
        //tracks in ecal
        List<Track> tracksInEcal = new ArrayList<Track>();
        //truth tracks
        List<Track> truthTracks = new ArrayList<Track>();

        //Map for storing Truth Tracks
        Map<Track,MCParticle> truthTracksMap = new HashMap<Track,MCParticle>();
        //List of Tracks that have been truth matched to a cluster
        List<Track> truthTracks_w_truthClusters = new ArrayList<Track>();
        //Map for storing truth tracks with truth cluster pairs
        Map<Track, Cluster> truthTracktruthClusterMap = new HashMap<Track, Cluster>();

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

        //Looping Over all Tracks
        int eventNumber = event.getEventNumber();
        for(Track track : tracks){

            double trackT;
            if (this.trackCollectionName.contains("GBLTracks")){
                trackT = TrackUtils.getTrackTime(track, hitToStrips, hitToRotated);
            }
            else {
                TrackData trackdata = (TrackData) TrktoData.from(track);
                trackT = trackdata.getTrackTime();
            }
            //For GBL Tracks, there is a |10|ns track time cut applied in
            //reconstruciton, that is not applied for Kalman Tracks. I am enforcing
            //this cut here, for a fair comparison of GBL v KF
            if(Math.abs(trackT) > 10.0)
                continue;

            //check track quality
            int charge = -1* (int)Math.signum(track.getTrackStates().get(0).getOmega());
            double chi2 = track.getChi2();
            int ndf = track.getNDF();

            List<Double> trackpos = getTrackPositionAtEcal(track);
            double trackx = trackpos.get(0);
            double tracky = trackpos.get(1);
            double trackz = trackpos.get(2);

            //Check if track at Ecal is within Ecal acceptance
            //We only want to feed Tracks that are expected to leave a hit in
            //the Ecal to the matching algorithm.
            boolean inEcalAccept = isTrackInEcal(trackx,tracky);
            if(inEcalAccept == true)
                tracksInEcal.add(track);

            double[] trackP = TrackUtils.getTrackStateAtLocation(track,TrackState.AtLastHit).getMomentum();
            double trackPmag = Math.sqrt(Math.pow(trackP[0],2) + Math.pow(trackP[1],2) + Math.pow(trackP[2],2));
            
            double[] trackPfinal = track.getTrackStates().get(track.getTrackStates().size()-1).getMomentum();
            double[] trackPinitial = track.getTrackStates().get(0).getMomentum();
            double deltaPmag = Math.sqrt(Math.pow(trackPfinal[0]-trackPinitial[0],2) + Math.pow(trackPfinal[1]-trackPinitial[1],2) + Math.pow(trackPfinal[2]-trackPinitial[2],2));

            //If want to select only FEEs for analysis
            boolean selectFEEs = false;
            if(selectFEEs){
                if(trackPmag > 2.5 || trackPmag < 2.0)
                    continue;
            }

            if(charge < 0){
                plots1D.get(String.format("%s_ele_track_momentum",this.trackCollectionName)).fill(trackPmag);
                plots1D.get(String.format("%s_ele_delta_momentum",this.trackCollectionName)).fill(deltaPmag);
                plots1D.get(String.format("%s_ele_track_chi2divndf",this.trackCollectionName)).fill(chi2/((double) ndf));
                plots2D.get(String.format("%s_ele_track_xypos_at_ecal",this.trackCollectionName)).fill(trackx,tracky);
            }
            else{
                plots1D.get(String.format("%s_pos_track_momentum",this.trackCollectionName)).fill(trackPmag);
                plots1D.get(String.format("%s_pos_delta_momentum",this.trackCollectionName)).fill(deltaPmag);
                plots1D.get(String.format("%s_pos_track_chi2divndf",this.trackCollectionName)).fill(chi2/((double) ndf));
                plots2D.get(String.format("%s_pos_track_xypos_at_ecal",this.trackCollectionName)).fill(trackx,tracky);
            }

            //TRUTH MATCH TRACKS TO MCPARTICLE
            MCParticle trackMCP =  newgetTrackMCP(event,track,this.trackCollectionName, true);
            if(trackMCP == null)
                continue;
            SimTrackerHit scoringplanehit = getTrackScoringPlaneHit(event, track, trackMCP, ecalScoringPlaneHitsCollectionName);
            if(scoringplanehit != null)
                trackScoringPlanePlots(event, track, scoringplanehit);

            //ONLY TRACKS THAT ARE TRUTH MATCHED TO MCPARTICLE PERSIST
            //BEYOND THIS POINT
            truthTracks.add(track);
            truthTracksMap.put(track, trackMCP);
            double mcpPmag = trackMCP.getMomentum().magnitude(); 
            double trackMCPratio = trackPmag/mcpPmag;

            if(charge < 0){
                plots1D.get(String.format("%s_ele_truth_track_momentum",this.trackCollectionName)).fill(trackPmag);
                plots1D.get(String.format("%s_ele_truth_track_MCP_momentum_ratio",this.trackCollectionName)).fill(trackMCPratio);
                plots2D.get(String.format("%s_ele_truth_track_xypos_at_ecal",this.trackCollectionName)).fill(trackx,tracky);
                plots2D.get(String.format("%s_ele_mcp_v_truth_track_momentum",this.trackCollectionName)).fill(mcpPmag,trackPmag);
            }
            else{
                plots1D.get(String.format("%s_pos_truth_track_momentum",this.trackCollectionName)).fill(trackPmag);
                plots1D.get(String.format("%s_pos_truth_track_MCP_momentum_ratio",this.trackCollectionName)).fill(trackMCPratio);
                plots2D.get(String.format("%s_pos_truth_track_xypos_at_ecal",this.trackCollectionName)).fill(trackx,tracky);
                plots2D.get(String.format("%s_pos_mcp_v_truth_track_momentum",this.trackCollectionName)).fill(trackPmag,mcpPmag);
            }

            //TRUTH MATCH TRACK TO CLUSTER BY MATCHING MCPARTICLES
            Cluster truthTrackCluster = null;
            List<Cluster> matchedClusters = new ArrayList<Cluster>();
            for(Map.Entry<Cluster,MCParticle> entry : truthClustersMap.entrySet()){
                MCParticle clusterMCP = entry.getValue();
                if(clusterMCP == trackMCP) {
                    truthTrackCluster = entry.getKey();
                    break;
                }
            }
            /*
            for(Map.Entry<Cluster,MCParticle> entry : truthClustersMap.entrySet()){
                MCParticle clusterMCP = entry.getValue();
                if(clusterMCP == trackMCP) {
                    matchedClusters.add(entry.getKey());
                    System.out.println("Track with momentum " + trackPmag + " matched to cluster with Energy " + entry.getKey().getEnergy());
                    System.out.println("This Cluster MCP has energy: " + clusterMCP.getEnergy());
                    System.out.println("This Track MCP has energy: " + trackMCP.getEnergy());
                    //truthTrackCluster = entry.getKey();                    
                }
            }

            Cluster bestMatch = null;
            double smallest_dr = 99999;
            if(matchedClusters == null || matchedClusters.size() < 1)
                continue;
            if(matchedClusters.size() > 1){
                System.out.println("Track has been truth matched to mutiple Clusters...");
                for(Cluster cluster : matchedClusters){
                    double clusterx = cluster.getPosition()[0];
                    double clustery = cluster.getPosition()[1];
                    double clusterz = cluster.getPosition()[2];
                    double clusterEnergy = cluster.getEnergy();
                    double dx = clusterx - trackx;
                    double dy = clustery - tracky;
                    double dz = clusterz - trackz;
                    double dr = Math.sqrt(Math.pow(dx,2) + Math.pow(dy,2) + Math.pow(dz,2));
                    double clustertime = ClusterUtilities.getSeedHitTime(cluster);
                    double dt = clustertime - trackClusterTimeOffset - trackT + 4.0;

                    //Defines Track top or bottom
                    double tanlambda = track.getTrackParameter(4);

                    //Check that Top/Bottom Tracks is only matched to
                    //Top/Bottom Clusters
                    if(clustery > 0 && tanlambda < 0)
                        continue;
                    if(clustery < 0 && tanlambda > 0)
                        continue;
                    if(Math.abs(dt) > 15)
                        continue;


                    //if(Math.abs(dt) > 4 || Math.abs(dx) > 20 || Math.abs(dy) > 15 || clusterEnergy/trackPmag < 0.6)
                       // continue;

                    plots1D.get(String.format("testing_track_w_multiple_cluster_matches_dx",trackCollectionName)).fill(clusterx - trackx);
                    plots1D.get(String.format("testing_track_w_multiple_cluster_matches_dy",trackCollectionName)).fill(clustery - tracky);
                    plots1D.get(String.format("testing_track_w_multiple_cluster_matches_dz",trackCollectionName)).fill(clusterz - trackz);
                    plots1D.get(String.format("testing_track_w_multiple_cluster_matches_dr",trackCollectionName)).fill(dr);
                    plots1D.get(String.format("testing_track_w_multiple_cluster_matches_dt",trackCollectionName)).fill(dt);
                    plots1D.get(String.format("testing_track_w_multiple_cluster_matches_EdivP",trackCollectionName)).fill(clusterEnergy/trackPmag);
                    plots1D.get(String.format("testing_track_w_multiple_cluster_matches_clusterMCP_Eratio",this.trackCollectionName)).fill(clusterEnergy/trackMCP.getEnergy());

                    if(dr < smallest_dr){
                        smallest_dr = dr;
                        bestMatch = cluster;
                    }
                    
                }
            }
            */

            //ONLY TRACKS THAT HAVE BEEN TRUTH MATCHED TO A CLUSTER PERSIST
            //BEYOND THIS POINT
            if(truthTrackCluster == null)
                continue;

            //Defines Track top or bottom
            double tanlambda = track.getTrackParameter(4);
            //Check that Top/Bottom Tracks is only matched to
            //Top/Bottom Clusters
            //if(truthTrackCluster.getPosition()[1] > 0 && tanlambda < 0)
            //    continue;
            //if(truthTrackCluster.getPosition()[1] < 0 && tanlambda > 0)
            //    continue;

            truthTracks_w_truthClusters.add(track);
            truthTracktruthClusterMap.put(track, truthTrackCluster);

            double clusterEnergy = truthTrackCluster.getEnergy();
            double clusterx = truthTrackCluster.getPosition()[0];
            double clustery = truthTrackCluster.getPosition()[1];
            double clusterz = truthTrackCluster.getPosition()[2];
            double truthTrackPmag = trackMCP.getMomentum().magnitude(); 

            //LOOKING AT ALL TRUTH TRACK CLUSTER PAIRS
            if(charge < 0){
                NtruthEleClustPairs = NtruthEleClustPairs + 1;
                plots1D.get(String.format("%s_ele_truth_track_cluster_pair_momentum",this.trackCollectionName)).fill(trackPmag);
                plots1D.get(String.format("%s_ele_truth_track_cluster_pair_px",this.trackCollectionName)).fill(trackP[1]);
                plots1D.get(String.format("%s_ele_truth_track_cluster_pair_py",this.trackCollectionName)).fill(trackP[2]);
                plots1D.get(String.format("%s_ele_truth_track_cluster_pair_pz",this.trackCollectionName)).fill(trackP[0]);
                plots1D.get(String.format("%s_ele_truth_track_cluster_pair_dx",trackCollectionName)).fill(clusterx - trackx);
                plots1D.get(String.format("%s_ele_truth_track_cluster_pair_dy",trackCollectionName)).fill(clustery - tracky);
                plots1D.get(String.format("%s_ele_truth_track_cluster_pair_dz",trackCollectionName)).fill(clusterz - trackz);
                plots1D.get(String.format("%s_ele_truth_track_cluster_pair_EdivP",trackCollectionName)).fill(clusterEnergy/trackPmag);
                plots1D.get(String.format("%s_ele_truth_track_cluster_pair_clusterMCP_Eratio",this.trackCollectionName)).fill(clusterEnergy/trackMCP.getEnergy());

                plots2D.get(String.format("%s_ele_truth_track_cluster_pair_xypos",this.trackCollectionName)).fill(trackx,tracky);
                
            }

            if(charge > 0){
                NtruthPosClustPairs = NtruthPosClustPairs + 1;
                plots1D.get(String.format("%s_pos_truth_track_cluster_pair_momentum",this.trackCollectionName)).fill(trackPmag);
                plots1D.get(String.format("%s_pos_truth_track_cluster_pair_px",this.trackCollectionName)).fill(trackP[1]);
                plots1D.get(String.format("%s_pos_truth_track_cluster_pair_py",this.trackCollectionName)).fill(trackP[2]);
                plots1D.get(String.format("%s_pos_truth_track_cluster_pair_pz",this.trackCollectionName)).fill(trackP[0]);
                plots1D.get(String.format("%s_pos_truth_track_cluster_pair_dx",trackCollectionName)).fill(clusterx - trackx);
                plots1D.get(String.format("%s_pos_truth_track_cluster_pair_dy",trackCollectionName)).fill(clustery - tracky);
                plots1D.get(String.format("%s_pos_truth_track_cluster_pair_EdivP",trackCollectionName)).fill(clusterEnergy/trackPmag);
                plots1D.get(String.format("%s_pos_truth_track_cluster_pair_dz",trackCollectionName)).fill(clusterz - trackz);
                plots1D.get(String.format("%s_pos_truth_track_cluster_pair_clusterMCP_Eratio",this.trackCollectionName)).fill(clusterEnergy/trackMCP.getEnergy());

                plots2D.get(String.format("%s_pos_truth_track_cluster_pair_xypos",this.trackCollectionName)).fill(trackx,tracky);
            }

            //LOOKING AT TRUTH TRACK CLUSTER PAIRS IN ECAL ACCEPTANCE
            if(inEcalAccept == true){
                if(charge < 0){
                    plots1D.get(String.format("%s_ele_truth_track_cluster_pair_insideEcal_momentum",this.trackCollectionName)).fill(trackPmag);
                    plots1D.get(String.format("%s_ele_truth_track_cluster_pair_insideEcal_dx",trackCollectionName)).fill(clusterx - trackx);
                    plots1D.get(String.format("%s_ele_truth_track_cluster_pair_insideEcal_dy",trackCollectionName)).fill(clustery - tracky);
                    plots1D.get(String.format("%s_ele_truth_track_cluster_pair_insideEcal_dz",trackCollectionName)).fill(clusterz - trackz);
                    plots1D.get(String.format("%s_ele_truth_track_cluster_pair_insideEcal_EdivP",trackCollectionName)).fill(clusterEnergy/trackPmag);
                    plots1D.get(String.format("%s_ele_truth_track_cluster_pair_insideEcal_clusterMCP_Eratio",this.trackCollectionName)).fill(clusterEnergy/trackMCP.getEnergy());

                    for(int i =0; i < eventNumber; i++){
                        plots2D.get(String.format("%s_ele_truth_track_cluster_pair_insideEcal_xypos",this.trackCollectionName)).fill(trackx,tracky);
                        plots2D.get(String.format("%s_ele_truth_track_cluster_pair_insideEcal_cluster_xypos",this.trackCollectionName)).fill(clusterx,clustery);
                    }
                    
                }

                if(charge > 0){
                    plots1D.get(String.format("%s_pos_truth_track_cluster_pair_insideEcal_momentum",this.trackCollectionName)).fill(trackPmag);
                    plots1D.get(String.format("%s_pos_truth_track_cluster_pair_insideEcal_dx",trackCollectionName)).fill(clusterx - trackx);
                    plots1D.get(String.format("%s_pos_truth_track_cluster_pair_insideEcal_dy",trackCollectionName)).fill(clustery - tracky);
                    plots1D.get(String.format("%s_pos_truth_track_cluster_pair_insideEcal_EdivP",trackCollectionName)).fill(clusterEnergy/trackPmag);
                    plots1D.get(String.format("%s_pos_truth_track_cluster_pair_insideEcal_dz",trackCollectionName)).fill(clusterz - trackz);
                    plots1D.get(String.format("%s_pos_truth_track_cluster_pair_insideEcal_clusterMCP_Eratio",this.trackCollectionName)).fill(clusterEnergy/trackMCP.getEnergy());

                    for(int i =0; i < eventNumber; i++){
                        plots2D.get(String.format("%s_pos_truth_track_cluster_pair_insideEcal_xypos",this.trackCollectionName)).fill(trackx,tracky);
                        plots2D.get(String.format("%s_pos_truth_track_cluster_pair_insideEcal_cluster_xypos",this.trackCollectionName)).fill(clusterx,clustery);
                    }
                
                }
            }



            //LOOKING AT TRACKS THAT ARE MATCHED TO CLUSTERS, BUT
            //WHERE TRACK IS NOT IN ECAL ACCEPTANCE
            if(inEcalAccept == false){
                int trackClusterTag = eventNumber;
                if(charge < 0){
                    plots1D.get(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_momentum",this.trackCollectionName)).fill(trackPmag);
                    plots1D.get(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_Pz",this.trackCollectionName)).fill(trackP[2]);
                    plots1D.get(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_dx",trackCollectionName)).fill(clusterx - trackx);
                    plots1D.get(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_dy",trackCollectionName)).fill(clustery - tracky);
                    plots1D.get(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_dz",trackCollectionName)).fill(clusterz - trackz);
                    plots1D.get(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_EdivP",trackCollectionName)).fill(clusterEnergy/trackPmag);
                    plots1D.get(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_clusterMCP_Eratio",this.trackCollectionName)).fill(clusterEnergy/trackMCP.getEnergy());

                    for(int i =0; i < trackClusterTag; i++){
                        plots2D.get(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_xypos",this.trackCollectionName)).fill(trackx,tracky);
                        plots2D.get(String.format("%s_ele_truth_track_cluster_pair_outsideEcal_cluster_xypos",this.trackCollectionName)).fill(clusterx,clustery);
                    }
                    
                }

                if(charge > 0){
                    plots1D.get(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_momentum",this.trackCollectionName)).fill(trackPmag);
                    plots1D.get(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_Pz",this.trackCollectionName)).fill(trackP[2]);
                    plots1D.get(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_dx",trackCollectionName)).fill(clusterx - trackx);
                    plots1D.get(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_dy",trackCollectionName)).fill(clustery - tracky);
                    plots1D.get(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_EdivP",trackCollectionName)).fill(clusterEnergy/trackPmag);
                    plots1D.get(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_dz",trackCollectionName)).fill(clusterz - trackz);
                    plots1D.get(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_clusterMCP_Eratio",this.trackCollectionName)).fill(clusterEnergy/trackMCP.getEnergy());

                    for(int i =0; i < trackClusterTag; i++){
                        plots2D.get(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_xypos",this.trackCollectionName)).fill(trackx,tracky);
                        plots2D.get(String.format("%s_pos_truth_track_cluster_pair_outsideEcal_cluster_xypos",this.trackCollectionName)).fill(clusterx,clustery);
                    }
                
                }
            }
        }


        //DEFINE VARIOUS METRICS TO EVALUATE PERFORMANCE
        for(int i = 0; i < truthTracks_w_truthClusters.size(); i++){
            int charge = -1* (int)Math.signum(truthTracks_w_truthClusters.get(i).getTrackStates().get(0).getOmega());
            if(charge < 0)
                ntruthpairs_ele = ntruthpairs_ele + 1;
            if(charge > 0)
                ntruthpairs_pos = ntruthpairs_pos + 1;
        }
        for(int i = 0; i < truthTracks.size(); i++){
            int charge = -1* (int)Math.signum(truthTracks.get(i).getTrackStates().get(0).getOmega());
            if(charge < 0)
                ntruthTracks_ele = ntruthTracks_ele + 1;
            if(charge > 0)
                ntruthTracks_pos = ntruthTracks_pos + 1;
        }
        for(int i = 0; i < tracks.size(); i++){
            int charge = -1* (int)Math.signum(tracks.get(i).getTrackStates().get(0).getOmega());
            if(charge < 0)
                nrecoTracks_ele = nrecoTracks_ele + 1;
            if(charge > 0)
                nrecoTracks_pos = nrecoTracks_pos + 1;
        }
        ntruthClusters = ntruthClusters + truthClusters.size();
        nClusters = nClusters + clusters.size();
        
        //Use matching algorithm to match tracks to clusters (without truth
        //info)
        Map<Track, Cluster> matchedTrackClusterMap = new HashMap<Track,Cluster>();
        //To check the performance of the Track+Cluster matching algorithm,
        //feed the algorithm the set of truth matched Tracks and truth
        //matched Clusters, so that comparison can be done.
        List<List<Track>> trackCollections = new ArrayList<List<Track>>();
        //trackCollections.add(tracksInEcal);
        trackCollections.add(tracks);
        matchedTrackClusterMap = matcher.matchTracksToClusters(event, trackCollections, clusters, cuts, -1, false, true, ecal, beamEnergy);
        if(matchedTrackClusterMap != null){

            for (Map.Entry<Track,Cluster> entry : matchedTrackClusterMap.entrySet()){

                Track track = entry.getKey();
                Cluster cluster = entry.getValue();

                int charge = -1* (int)Math.signum(track.getTrackStates().get(0).getOmega());
                double[] trackP = TrackUtils.getTrackStateAtLocation(track,TrackState.AtLastHit).getMomentum();
                double trackPmag = Math.sqrt(Math.pow(trackP[0],2) + Math.pow(trackP[1],2) + Math.pow(trackP[2],2));

                //plot xy positions
                List<Double> trackpos = getTrackPositionAtEcal(track);
                boolean isinEcal = isTrackInEcal(trackpos.get(0),trackpos.get(1)); //if track is outside of ecal, returns false

                boolean truthTrack = false;
                double truthTrackPDGID = -9999;
                MCParticle truthTrackMCP = null;
                boolean trackOutsideEcal = false;

                boolean truthCluster = false;
                double truthClusterPDGID = 9999;
                MCParticle truthClusterMCP = null;

                boolean noCluster = false;
                boolean fakeMatch = false;
                boolean unknownMatch = false;
                boolean goodMatch = false;
                boolean missedMatch = false;
                boolean noMatch = false;
                boolean noMatchExpected = false;


                if(cluster == null)
                    noCluster = true;

                if(!isinEcal)
                    trackOutsideEcal = true;

                if(noCluster == false){
                    if(truthTracksMap.containsKey(track)){
                        truthTrack = true;
                        truthTrackPDGID = truthTracksMap.get(track).getPDGID();
                        truthTrackMCP = truthTracksMap.get(track);
                    }
                    if(truthClustersMap.containsKey(cluster)){
                        truthCluster = true;
                        truthClusterPDGID = truthClustersMap.get(cluster).getPDGID();
                        truthClusterMCP = truthClustersMap.get(cluster);
                    }

                    if(truthTrack & truthCluster){
                        if(truthTrackMCP == truthClusterMCP)
                            goodMatch = true;
                        else
                            fakeMatch = true;
                    }

                    if(truthTrack & !truthCluster){
                        //if MCP has track and cluster, but track is paired to
                        //different cluster, bad match
                        if(truthClustersMap.containsValue(truthTrackMCP))
                            fakeMatch = true;
                        else
                            unknownMatch = true;
                    }

                    if(truthCluster & !truthTrack){
                        if(truthTracksMap.containsValue(truthClusterMCP))
                            fakeMatch = true;
                        else if(truthClusterPDGID == 22 || (truthClusterPDGID == -11 & charge < 0) || (truthClusterPDGID == 11 & charge > 0))                    {
                            fakeMatch = true;
                        }
                        else
                            unknownMatch = true;
                    }

                    if(!truthCluster & !truthTrack){
                        unknownMatch = true;
                    }
                }

                else if(noCluster == true){
                    noMatch = true;
                    if(truthTracksMap.containsKey(track)){
                        if(truthClustersMap.containsValue(truthTrackMCP))
                            missedMatch = true;
                    }
                }

                //add numbers
                
                if(noMatch & truthTrack){
                    if(charge < 0)
                        noMatchTruthTrack_ele = noMatchTruthTrack_ele + 1;
                    else
                        noMatchTruthTrack_pos = noMatchTruthTrack_pos + 1;
                }

                if(unknownMatch & truthTrack){
                    if(charge < 0)
                        unknownMatchTruthTrack_ele = unknownMatchTruthTrack_ele + 1;
                    else
                        unknownMatchTruthTrack_pos = unknownMatchTruthTrack_pos + 1;
                }

                if(unknownMatch & truthCluster){
                    if(charge < 0)
                        unknownMatchTruthCluster_ele = unknownMatchTruthCluster_ele + 1;
                    else
                        unknownMatchTruthCluster_pos = unknownMatchTruthCluster_pos + 1;
                }

                if(truthTrack & trackOutsideEcal){
                    if(charge < 0)
                        truthTrackOutsideEcal_ele = truthTrackOutsideEcal_ele + 1;
                    else
                        truthTrackOutsideEcal_pos = truthTrackOutsideEcal_pos + 1;
                }

                if(trackOutsideEcal){
                    if(charge < 0)
                        trackOutsideEcal_ele = trackOutsideEcal_ele + 1;
                    else
                        trackOutsideEcal_pos = trackOutsideEcal_pos + 1;
                }

                if(goodMatch){
                    if(charge < 0)
                        goodMatch_ele = goodMatch_ele + 1;
                    else
                        goodMatch_pos = goodMatch_pos + 1;
                }

                if(fakeMatch){
                    if(charge < 0)
                        fakeMatch_ele = fakeMatch_ele + 1;
                    else
                        fakeMatch_pos = fakeMatch_pos + 1;
                }

                if(noMatch){
                    if(charge < 0)
                        noMatch_ele = noMatch_ele + 1;
                    else
                        noMatch_pos = noMatch_pos + 1;
                }

                if(unknownMatch){
                    if(charge < 0)
                        unknownMatch_ele = unknownMatch_ele + 1;
                    else
                        unknownMatch_pos = unknownMatch_pos + 1;
                }

                if(missedMatch){
                    if(charge < 0)
                        missedMatch_ele = missedMatch_ele + 1;
                    else
                        missedMatch_pos = missedMatch_pos + 1;
                }
            }
        }

        List<MCParticle> possibleTrackMCPs = getPossibleTrackMCPs(event, 9);
        NpossibleTracks = NpossibleTracks + possibleTrackMCPs.size();
        NrecoTruthTracks = NrecoTruthTracks + truthTracksMap.size();
        Map<MCParticle,int[]> nDuplicateMCPs = new HashMap<MCParticle, int[]>();

        for(Map.Entry<Track, MCParticle> entry : truthTracksMap.entrySet()){
            MCParticle particle = entry.getValue(); 
            if(!nDuplicateMCPs.containsKey(particle)){
                nDuplicateMCPs.put(particle, new int[1]);
                nDuplicateMCPs.get(particle)[0] = 0;
            }

            nDuplicateMCPs.get(particle)[0]++;
        }
        for(Map.Entry<MCParticle,int[]> entry : nDuplicateMCPs.entrySet()){
            int nentries = entry.getValue()[0];
            List<Track> dupTracks = new ArrayList<Track>();
            if(nentries > 1){
                for(Map.Entry<Track, MCParticle> entry2 : truthTracksMap.entrySet())
                    if(entry2.getValue() == entry.getKey())
                        dupTracks.add(entry2.getKey());
                for(Track track : dupTracks){
                    int charge = -1* (int)Math.signum(track.getTrackStates().get(0).getOmega());
                    double[] trackP = track.getMomentum();
                    double trackPmag = Math.sqrt(Math.pow(trackP[0],2) + Math.pow(trackP[1],2) + Math.pow(trackP[2],2));
                    if(charge < 0)
                        plots1D.get(String.format("%s_ele_duplicate_track_momentum",this.trackCollectionName)).fill(trackPmag);
                    else
                        plots1D.get(String.format("%s_pos_duplicate_track_momentum",this.trackCollectionName)).fill(trackPmag);
                }
                nDuplicates = nDuplicates + 1.0;
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

    public boolean isTrackInSlot(double trackx, double tracky){

        //Define first order Ecal geometry --> assumed square here,
        //smaller than actual beamgap. Improve x geometry 
        double ecalx1 = -276.0; //mm
        double ecalx2 = 361.0;
        double ecaly1 = 91.0;
        double ecaly2 = -91.0;
        double bgapup = 22.4;
        double bgapdown = -22.3;

        double eholex11 = -93.0;
        double eholex12 = -85.0;

        double eholex22 = 16.0;
        double eholex21 = 29;

        double eholey12 = 37.0;
        double eholey11 = 22.4;

        double eholey22 = -37.0;
        double eholey21 = -22.3;

        boolean inSlot = false;

        if((trackx > ecalx1 && trackx < ecalx2) && (tracky < ecaly1 && tracky > ecaly2) && tracky < bgapup && tracky > bgapdown)
            inSlot = true;
        return inSlot;
    }


    public SimTrackerHit getTrackScoringPlaneHit(EventHeader event, Track track, MCParticle trackMCP, String ecalScoringPlaneHitsCollectionName) {

        List<SimTrackerHit> scoringPlaneHits = event.get(SimTrackerHit.class, ecalScoringPlaneHitsCollectionName);


        //Check for simtrackerhit MCP that matches trackMCP
        if(trackMCP == null)
            return null;
        SimTrackerHit matchedScoringPlaneHit = null;
        for(SimTrackerHit scoringPlaneHit : scoringPlaneHits){
            // If the MC particles don't match, move on to the next particle
            if(!(scoringPlaneHit.getMCParticle() == trackMCP)) continue;
            matchedScoringPlaneHit = scoringPlaneHit;
            // Once a match is found, there is no need to loop through the rest of the list
            break;
        }
        return matchedScoringPlaneHit;
    }

    public void getTruthClustersMap(EventHeader event, List<Cluster> clusters){

        HashMap<MCParticle, HashSet<Cluster>> mcpClusterMap = new HashMap<MCParticle, HashSet<Cluster>>();
        for(Cluster cluster : clusters){
            double clustx = cluster.getPosition()[0];
            double clusty = cluster.getPosition()[1];
            double clustz = cluster.getPosition()[2];
            double clusterEnergy = cluster.getEnergy();

            MCParticle clusterMCP = getClusterMCP(cluster,event, false);
            if(clusterMCP == null)
                continue;
            //check status of clusterMCP
            System.out.println("clusterMCP PDGID = " + clusterMCP.getPDGID());
            if(clusterMCP.getSimulatorStatus().isBackscatter()){
                System.out.println("ERROR: Cluster MCP is BackScatter!");
                System.out.println("MCP origin z: " + clusterMCP.getOriginZ());
                plots2D.get("cluster_mcp_isBackscatter_zy").fill(clusterMCP.getOriginZ(), clusterMCP.getOriginY());
            }
            if(clusterMCP.getSimulatorStatus().isDecayedInTracker())
                System.out.println("Status Check: Cluster MCP is decayed in Tracker");
            if(clusterMCP.getSimulatorStatus().isDecayedInCalorimeter())
                System.out.println("Status Check: Cluster MCP is decayed in Calorimeter");
            if(clusterMCP.getSimulatorStatus().isCreatedInSimulation())
                System.out.println("Status Check: Cluster MCP is created in Simulation");
            else
                System.out.println("Status Check: Cluster MCP is NOT created in Simulation");

            HashSet<Cluster> c = new HashSet<Cluster>();
            if(mcpClusterMap.containsKey(clusterMCP)){
                c = mcpClusterMap.get(clusterMCP);
                c.add(cluster);
                mcpClusterMap.put(clusterMCP, c);
            }
            else{
                c.add(cluster);
                mcpClusterMap.put(clusterMCP, c);
            }
        }   

        for(Map.Entry<MCParticle, HashSet<Cluster>> entry : mcpClusterMap.entrySet()){
            MCParticle mcp = entry.getKey();
            HashSet<Cluster> tClusters = entry.getValue();
            System.out.println("Checking For Duplicates");
            if(tClusters.size() > 0)
                plots1D.get("nclusters_matched_to_single_mcp").fill(tClusters.size());

            if(tClusters.size() > 1){
                System.out.println("Multiple clusters matched to same MCP!");
                System.out.println(tClusters.size() + " different Clusters matched to this MCP");
                //check status of clusterMCP
                System.out.println("mcp PDGID = " + mcp.getPDGID());
                System.out.println("mcp OriginZ: " + mcp.getOriginZ());
                System.out.println("mcp EndpointX: " + mcp.getEndPoint().x());
                System.out.println("mcp EndpointY: " + mcp.getEndPoint().y());
                System.out.println("mcp EndpointZ: " + mcp.getEndPoint().z());

                if(mcp.getPDGID() == 22){
                    System.out.println("Daughters of duplicate photon cluster: ");
                    List<MCParticle> daughters = mcp.getDaughters();
                    for(MCParticle daughter : daughters){
                        System.out.println("Daughter PDGID: "+ daughter.getPDGID());
                        System.out.println("Daughter Energy: " + daughter.getEnergy());
                    }
                }

                if(mcp.getSimulatorStatus().isBackscatter()){
                    System.out.println("ERROR: Cluster MCP is BackScatter!");
                    System.out.println("MCP origin z: " + mcp.getOriginZ());
                }
                if(mcp.getSimulatorStatus().isDecayedInTracker())
                    System.out.println("Status Check: Cluster MCP is decayed in Tracker");
                if(mcp.getSimulatorStatus().isDecayedInCalorimeter())
                    System.out.println("Status Check: Cluster MCP is decayed in Calorimeter");
                if(mcp.getSimulatorStatus().isCreatedInSimulation())
                    System.out.println("Status Check: Cluster MCP is created in Simulation");
                else
                    System.out.println("Status Check: Cluster MCP is NOT created in Simulation");

                System.out.println("Duplicate MCP Energy: " + mcp.getEnergy());
                for(Cluster cluster : tClusters){
                    System.out.println("Duplicate Cluster Energy = " + cluster.getEnergy());
                    System.out.println("Duplicate Cluster ypos = " + cluster.getPosition()[1]);
                    System.out.println("Duplicate Cluster xpos = " + cluster.getPosition()[0]);

                    //MCParticle trash = getClusterMCP(cluster, event, true);

                    if(mcp.getPDGID() == 22){
                        plots2D.get("cluster_photon_mcp_duplicates_xy").fill(cluster.getPosition()[0],cluster.getPosition()[1]);
                    }
                    if(mcp.getPDGID() == 11)
                        plots2D.get("cluster_electron_mcp_duplicates_xy").fill(cluster.getPosition()[0],cluster.getPosition()[1]);
                    if(mcp.getPDGID() == -11)
                        plots2D.get("cluster_positron_mcp_duplicates_xy").fill(cluster.getPosition()[0],cluster.getPosition()[1]);
                }

            }

        }

    }
    

    public Map<Cluster,MCParticle> getUniqueTruthClusters(EventHeader event, List<Cluster> clusters){

        System.out.println("CHECK EVENT NUMBER: " + event.getEventNumber());

        getTruthClustersMap(event, clusters);
        if(clusters != null)
            plots1D.get("cluster_event_multiplicity").fill(clusters.size());
        else
            plots1D.get("cluster_event_multiplicity").fill(0);

        //Truth Match Clusters in LCIO to MCParticles    
        List<Cluster> truthClustersUnsorted = new ArrayList<Cluster>();
        Map<Cluster,MCParticle> truthClustersMapUnsorted = new HashMap<Cluster,MCParticle>();
        for(Cluster cluster : clusters){
            double clustx = cluster.getPosition()[0];
            double clusty = cluster.getPosition()[1];
            double clustz = cluster.getPosition()[2];
            double clusterEnergy = cluster.getEnergy();
            boolean edge = false;
            boolean notedge = false;

            //if(clustx > -280 && clustx < -267)
              //  edge = true;
            //if(clustx < 350 && clustx >337)
              //  edge = true;
            if((clustx < -85) & ((clusty > 0 & clusty < 30) || (clusty < 0 & clusty > -30 ) ))
                edge = true;
            if((clustx > 15) & ((clusty > 0 & clusty < 30) || (clusty < 0 & clusty > -30 )))
                edge = true;
            if((clustx > -85 & clustx < 15) & ((clusty > 0 & clusty < 43) || (clusty < 0 & clusty > -44 )))
                edge = true;

            if((clustx > 268 & clustx < -85) & ((clusty > 32) || (clusty < -32 ) ))
                notedge = true;
            if((clustx > 15 & clustx < 353) & ((clusty > 32) || (clusty < -32 )))
                notedge = true;
            if((clustx > -85 & clustx < 15) & ((clusty > 45) || (clusty < -46 )))
                notedge = true;
            

            MCParticle clusterMCP = getClusterMCP(cluster,event,true);
            plots2D.get(String.format("ecal_cluster_positions_xy_plane")).fill(clustx,clusty);
            plots1D.get(String.format("cluster_energy",trackCollectionName)).fill(clusterEnergy);
            if(clusterMCP == null)
                continue;

            if(clusterMCP.getPDGID() == 22){
                plots1D.get("cluster_photon_mcp_truth_originz").fill(clusterMCP.getOriginZ());
                plots1D.get("cluster_photon_mcp_truth_energy").fill(clusterMCP.getEnergy());
                plots2D.get(String.format("cluster_vs_photon_mcp_truth_energy")).fill(clusterEnergy,clusterMCP.getEnergy());
                if(edge)
                    plots2D.get(String.format("cluster_vs_photon_mcp_truth_energy_onEdge")).fill(clusterEnergy,clusterMCP.getEnergy());
                if(notedge)
                    plots2D.get(String.format("cluster_vs_photon_mcp_truth_energy_offEdge")).fill(clusterEnergy,clusterMCP.getEnergy());
                else if(!edge & !notedge)
                    plots2D.get(String.format("cluster_vs_photon_mcp_truth_energy_unsureEdge")).fill(clusterEnergy,clusterMCP.getEnergy());
            }
            if(clusterMCP.getPDGID() == 11){
                plots1D.get("cluster_electron_mcp_truth_originz").fill(clusterMCP.getOriginZ());
                plots1D.get("cluster_electron_mcp_truth_energy").fill(clusterMCP.getEnergy());
                plots2D.get(String.format("cluster_vs_electron_mcp_truth_energy")).fill(clusterEnergy,clusterMCP.getEnergy());
                if(edge)
                    plots2D.get(String.format("cluster_vs_electron_mcp_truth_energy_onEdge")).fill(clusterEnergy,clusterMCP.getEnergy());
                if(notedge)
                    plots2D.get(String.format("cluster_vs_electron_mcp_truth_energy_offEdge")).fill(clusterEnergy,clusterMCP.getEnergy());
                else if(!edge & !notedge)
                    plots2D.get(String.format("cluster_vs_electron_mcp_truth_energy_unsureEdge")).fill(clusterEnergy,clusterMCP.getEnergy());

            }
            if(clusterMCP.getPDGID() == -11){
                plots1D.get("cluster_positron_mcp_truth_originz").fill(clusterMCP.getOriginZ());
                plots1D.get("cluster_positron_mcp_truth_energy").fill(clusterMCP.getEnergy());
                plots2D.get(String.format("cluster_vs_positron_mcp_truth_energy")).fill(clusterEnergy,clusterMCP.getEnergy());
                if(edge)
                    plots2D.get(String.format("cluster_vs_positron_mcp_truth_energy_onEdge")).fill(clusterEnergy,clusterMCP.getEnergy());
                if(notedge)
                    plots2D.get(String.format("cluster_vs_positron_mcp_truth_energy_offEdge")).fill(clusterEnergy,clusterMCP.getEnergy());
                else if(!edge & !notedge)
                    plots2D.get(String.format("cluster_vs_positron_mcp_truth_energy_unsureEdge")).fill(clusterEnergy,clusterMCP.getEnergy());

            }

            double ecalFacePosZ = 1330;
            if(clusterMCP.getEndPoint().y() >= ecalFacePosZ){
                plots2D.get(String.format("clustery_v_mcpy")).fill(clusty,clusterMCP.getEndPoint().y());
                plots2D.get(String.format("clusterx_v_mcpx")).fill(clustx,clusterMCP.getEndPoint().x());
            }

            truthClustersUnsorted.add(cluster);
            truthClustersMapUnsorted.put(cluster, clusterMCP);

        }
        plots1D.get("cluster_mcp_truth_match_event_multiplicity").fill(truthClustersUnsorted.size());

        //Check distributions of MCP matched to multiple clusters
        List<Cluster> skipentry = new ArrayList<Cluster>();
        for(Map.Entry<Cluster, MCParticle> c : truthClustersMapUnsorted.entrySet()){
            Cluster cluster1 = c.getKey();
            if(skipentry.contains(cluster1))
                continue;
            skipentry.add(cluster1);
            MCParticle mcp1 = c.getValue();
            double mcp1Energy = mcp1.getEnergy();
            double cluster1Energy = cluster1.getEnergy();
            for(Map.Entry<Cluster, MCParticle> cc : truthClustersMapUnsorted.entrySet()){
                Cluster cluster2 = cc.getKey();
                if(skipentry.contains(cluster2))
                    continue;
                MCParticle mcp2 = cc.getValue();
                double mcp2Energy = mcp2.getEnergy();
                double cluster2Energy = cluster2.getEnergy();
                double dr = Math.sqrt(Math.pow(cluster1.getPosition()[0]-cluster2.getPosition()[0],2) + Math.pow(cluster1.getPosition()[1]-cluster2.getPosition()[1],2));
                if(mcp1 == mcp2){
                    if(mcp1.getPDGID() == 22){
                        plots1D.get("cluster_photon_mcp_duplicates_dx").fill((cluster1.getPosition()[0]-cluster2.getPosition()[0]));
                        plots1D.get("cluster_photon_mcp_duplicates_dy").fill((cluster1.getPosition()[1]-cluster2.getPosition()[1]));
                        plots1D.get("cluster_photon_mcp_duplicates_dr").fill(dr);
                        plots1D.get("cluster_photon_mcp_duplicates_dE").fill(Math.abs(cluster1Energy-cluster2Energy));
                        plots2D.get("cluster_photon_mcp_duplicates_dx_v_dy").fill(cluster1.getPosition()[0]-cluster2.getPosition()[0],cluster1.getPosition()[1]-cluster2.getPosition()[1]);
                        plots1D.get("cluster_photon_mcp_duplicates_originz").fill(mcp1.getOriginZ());
                    }
                    if(mcp1.getPDGID() == 11){
                        plots1D.get("cluster_electron_mcp_duplicates_dx").fill((cluster1.getPosition()[0]-cluster2.getPosition()[0]));
                        plots1D.get("cluster_electron_mcp_duplicates_dy").fill((cluster1.getPosition()[1]-cluster2.getPosition()[1]));
                        plots1D.get("cluster_electron_mcp_duplicates_dE").fill(Math.abs(cluster1Energy-cluster2Energy));
                        plots1D.get("cluster_electron_mcp_duplicates_dr").fill(dr);
                        plots2D.get("cluster_electron_mcp_duplicates_dx_v_dy").fill(cluster1.getPosition()[0]-cluster2.getPosition()[0],cluster1.getPosition()[1]-cluster2.getPosition()[1]);
                        plots1D.get("cluster_electron_mcp_duplicates_originz").fill(mcp1.getOriginZ());
                    }
                    if(mcp1.getPDGID() == -11){
                        plots1D.get("cluster_positron_mcp_duplicates_dx").fill((cluster1.getPosition()[0]-cluster2.getPosition()[0]));
                        plots1D.get("cluster_positron_mcp_duplicates_dy").fill((cluster1.getPosition()[1]-cluster2.getPosition()[1]));
                        plots1D.get("cluster_positron_mcp_duplicates_dE").fill(Math.abs(cluster1Energy-cluster2Energy));
                        plots1D.get("cluster_positron_mcp_duplicates_dr").fill(dr);
                        plots2D.get("cluster_positron_mcp_duplicates_dx_v_dy").fill(cluster1.getPosition()[0]-cluster2.getPosition()[0],cluster1.getPosition()[1]-cluster2.getPosition()[1]);
                        plots1D.get("cluster_positron_mcp_duplicates_originz").fill(mcp1.getOriginZ());
                    }
                    
                }
            }
        }
            
        /*

        //Duplicate clusters are truth-matched to the same MCParticle.
        //To remove duplicates, loop over map and keep cluster_mcp pair
        //with the smallest energy difference
        Map<Cluster,MCParticle> truthClustersMap = new HashMap<Cluster,MCParticle>();
        List<Cluster> truthClusters = new ArrayList<Cluster>();

        List<Cluster> skipentry = new ArrayList<Cluster>();
        for(Map.Entry<Cluster, MCParticle> c : truthClustersMapUnsorted.entrySet()){
            if(skipentry.contains(c.getKey()))
                continue;
            skipentry.add(c.getKey());
            MCParticle mcp1 = c.getValue();
            double mcpenergy = mcp1.getEnergy();
            double cenergy = c.getKey().getEnergy();
            double dc = Math.abs(cenergy - mcpenergy);
            double smallestdE = dc;
            Cluster smallestdEcluster = c.getKey();
            for(Map.Entry<Cluster, MCParticle> cc : truthClustersMapUnsorted.entrySet()){
                if(skipentry.contains(cc.getKey()))
                    continue;
                if(mcp1 == cc.getValue()){
                    skipentry.add(cc.getKey());
                    double ccenergy = cc.getKey().getEnergy();
                    double dcc = Math.abs(ccenergy - mcpenergy);
                    if(dcc < smallestdE){
                        smallestdE = dcc;
                        smallestdEcluster = cc.getKey();
                    }
                }
            }

            truthClustersMap.put(smallestdEcluster,truthClustersMapUnsorted.get(smallestdEcluster));
            truthClusters.add(smallestdEcluster);
        }

        return truthClustersMap;
        */

        Map<Cluster,MCParticle> truthClustersMap = new HashMap<Cluster,MCParticle>();
        List<Cluster> truthClusters = new ArrayList<Cluster>();
        for(Map.Entry<Cluster, MCParticle> c : truthClustersMapUnsorted.entrySet()){
            boolean hasDup = false;
            Cluster cluster1 = c.getKey();
            MCParticle mcp1 = c.getValue();
            for(Map.Entry<Cluster, MCParticle> cc : truthClustersMapUnsorted.entrySet()){
                Cluster cluster2 = cc.getKey();
                MCParticle mcp2 = cc.getValue();
                if(cluster1 != cluster2 && mcp1 == mcp2)
                    hasDup = true;
            }
            if(hasDup == false)
                truthClustersMap.put(cluster1,mcp1);
        }
        return truthClustersMap;
    }


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
        double[] trackP = track.getTrackStates().get(track.getTrackStates().size()-1).getMomentum();
        double trackPmag = Math.sqrt(Math.pow(trackP[0],2) + Math.pow(trackP[1],2) + Math.pow(trackP[2],2));
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


    public List<Double> extrapTrackToCrystalCenter(Track track){

        List<Double> trackpos = getTrackPositionAtEcal(track);
        double trackx = trackpos.get(0);
        double tracky = trackpos.get(1);
        double trackz = trackpos.get(2);
        double[] trackP = track.getTrackStates().get(track.getTrackStates().size()-1).getMomentum();
        double px = trackP[0];
        double py = trackP[1];
        double pz = trackP[2];
        double crystalDepth = 80.0/1000;

        double deltax;
        double deltay;
        for(int i = 0; i < crystalDepth; i++){
            double iter = (double) i;
            deltax = (px/pz)*iter;
            deltay = (py/pz)*iter;
            trackx = trackx + deltax*1000;
            tracky = tracky + deltay*1000;
            trackz = trackz + iter;
        }

        List<Double> newTrackPos = new ArrayList<Double>();
        newTrackPos.add(trackx);
        newTrackPos.add(tracky);
        newTrackPos.add(trackz);

        return newTrackPos;
    }

/**
     * Get the MC particle associated with a track.
     * 
     * @param track : Track to get the MC particle for
     * @return The MC particle associated with the track
     */

    public void printMCPtree(EventHeader event){
        List<MCParticle> mcps = event.get(MCParticle.class,"MCParticle");

        for(int i = 0; i < mcps.size(); i++){
            MCParticle mcp = mcps.get(i);
            System.out.println("MCParticle " + i + " PDGID: " + mcp.getPDGID());
            System.out.println("MCParticle " + i + " Energy: " + mcp.getEnergy());
            System.out.println("MCParticle " + i + " Momentum: " + mcp.getPX() + "," + mcp.getPY() + "," + mcp.getPZ());
            System.out.println("MCParticle " + i + " Origin: " + mcp.getOriginX() + "," + mcp.getOriginY() + "," + mcp.getOriginZ());
            //if(mcp.getEndPoint())
            try{
                System.out.println("MCParticle " + i + " Endpoint: " + mcp.getEndPoint().x() + "," + mcp.getEndPoint().y() + "," + mcp.getEndPoint().z());
            }
            catch (Throwable t){
                System.out.println("No endpoint available");
            }
            if(mcp.getParents().size() != 0){   
                System.out.println("MCParticle " + i + " Parent PDGID: " + mcp.getParents().get(0).getPDGID());
                System.out.println("MCParticle " + i + " Parent Energy: " + mcp.getParents().get(0).getEnergy());
                System.out.println("MCParticle " + i + " Parent Momentum: " + mcp.getParents().get(0).getPX() + "," + mcp.getParents().get(0).getPY() + "," + mcp.getParents().get(0).getPZ());
            }
            else
                System.out.println("No parents of this MCP");
        }

    }

    private List<MCParticle> getMCParticlesFromLCIO(EventHeader event){
        List<MCParticle> mcparticles = event.get(MCParticle.class, "MCParticle");
        return mcparticles;
    }


    private Map<Cluster, MCParticle> getClusterMcpMap(List<Cluster> clusters, EventHeader event, boolean verbose){

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

                //Test if we should choose largest MCParticle in the simcalhit,
                //or parent of largest contributor 03021515
                /*
                if(largest_simcalhit.getContributedEnergy(i) > bestMCPEnergyContr){
                    bestMCPEnergyContr = largest_simcalhit.getContributedEnergy(i);
                    bestMCP = largest_simcalhit.getMCParticle(i);
                }
                */
            }
            /*
            //if MCP energy is less than what it contributes...it's likely that
            //the energy contribution assignments were wrong in simulation, and
            //the simcalhit should be matched to the parent of that contributor
            if(bestMCPEnergyContr > bestMCP.getEnergy()){
                if(verbose)
                    System.out.println("MCParticle with energy " + bestMCP.getEnergy() + " claims to have contributed " + bestMCPEnergyContr + "to simcalorimeter hit" );
                if(bestMCP.getParents() != null){
                    if(verbose)
                        System.out.println("MCParticle has parents in this simcalorimeter hit");
                    double bestParentEnergy = 0.0;
                    MCParticle bestParent = null;
                    for(MCParticle parent : bestMCP.getParents()){
                        double energy = parent.getEnergy();
                        if(verbose)
                            System.out.println("Parent energy: " + energy);
                        if(energy > bestParentEnergy){
                            bestParentEnergy = energy;
                            bestParent = parent;
                        }
                    }
                    if(verbose)
                        System.out.println("Matching parent with energy " + bestParentEnergy + "instead of previous match with energy " + bestMCP.getEnergy());

                    bestMCP = bestParent;
                }
                else{
                    if(verbose)
                        System.out.println("No parents found for this strange case....ignore this match");
                    bestMCP = null;
                }
            }
                */

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
            plots1D.get("nclusters_matched_to_single_mcp").fill(tClusters.size());

            if(tClusters.size() < 2){
                for(Cluster tCluster : tClusters){
                    mcpClusterMapFinal.put(tCluster, mcp);

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
                }

                if(recoveredDup.size() > 1){
                    for(int i =0; i < recoveredDup.size()-1; i++){
                        plots1D.get("cluster_truth_stage_3_cut_remaining_duplicate_mcp_match_dx").fill(recoveredDup.get(i).getPosition()[0] - recoveredDup.get(i+1).getPosition()[0]);
                        plots1D.get("cluster_truth_stage_3_cut_remaining_duplicate_mcp_match_dy").fill(recoveredDup.get(i).getPosition()[1] - recoveredDup.get(i+1).getPosition()[1]);

                    }
                }
                else if(recoveredDup.size() > 0 & recoveredDup.size() < 2) {
                    mcpClusterMapFinal.put(recoveredDup.get(0),mcp);
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
        return mcpClusterMapFinal;
    }

    private MCParticle getClusterMCP(Cluster cluster, EventHeader event, boolean verbose) {

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
            return null;
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
        double bestMCPContrE = 0.0;
        if(largest_simcalhit == null)
            return null;
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
                double distance = Math.sqrt( Math.pow(largest_simcalhit.getPosition()[0] - mcp.getEndPoint().x(),2) + Math.pow(largest_simcalhit.getPosition()[1] - mcp.getEndPoint().y(),2));
            }

            if(largest_simcalhit.getContributedEnergy(i) > bestMCPEnergyContr){
                bestMCPEnergyContr = largest_simcalhit.getContributedEnergy(i);
                bestMCP = largest_simcalhit.getMCParticle(i);
                bestMCPContrE = largest_simcalhit.getContributedEnergy(i);
            }
        }

        if(bestMCP == null)
            return null;

        double ecalFacePosZ = 1330.0; //mm
        if(bestMCP.getOriginZ() > ecalFacePosZ)
            return null;

        double distance = Math.sqrt( Math.pow(largest_simcalhit.getPosition()[0] - bestMCP.getEndPoint().x(),2) + Math.pow(largest_simcalhit.getPosition()[1] - bestMCP.getEndPoint().y(),2));

        double energyRatio = cluster.getEnergy()/bestMCP.getEnergy();
        plots1D.get(String.format("cluster_truthMCP_energy_ratio",this.trackCollectionName)).fill(energyRatio);

        return bestMCP;
    }

    private SimTrackerHit getBestSimHit(RawTrackerHit rawhit, RelationalTable rawtomc){

        if(rawhit == null)
            return null;
        Set<SimTrackerHit> simhits = rawtomc.allFrom(rawhit);
        SimTrackerHit rawhitSimhit = null;
        double simhitMaxE = 0.0;
        for(SimTrackerHit simhit : simhits){
            if (simhit != null && simhit.getMCParticle() != null) {
                double simhitEnergy = simhit.getdEdx();
                if(simhitEnergy > simhitMaxE){
                    simhitMaxE = simhitEnergy;
                    rawhitSimhit = simhit;
                }
            }
        }
        return rawhitSimhit;
     
    }

    private Set<RawTrackerHit> getRawHitsFromTrackerHit(TrackerHit hit, RelationalTable rawtomc){

        Map<RawTrackerHit, SimTrackerHit> largestHitOnLayerMap = new HashMap<RawTrackerHit, SimTrackerHit>();
        Map<RawTrackerHit, Integer> layerMap = new HashMap<RawTrackerHit, Integer>(); 
        Set<RawTrackerHit> rawhitsPerLayer = new HashSet<RawTrackerHit>();
        Set<Integer> layers = new HashSet<Integer>();
        for(RawTrackerHit rawhit : (List<RawTrackerHit>) hit.getRawHits()){
            int layer = rawhit.getLayerNumber(); 
            layerMap.put(rawhit, layer);
            layers.add(layer);
        }

        for(int layer : layers){
            List<RawTrackerHit> rawhitsonLayer = new ArrayList<RawTrackerHit>();
            for(Map.Entry<RawTrackerHit,Integer> entry : layerMap.entrySet()){
                if(entry.getValue().equals(layer)){
                    rawhitsonLayer.add(entry.getKey());
                }
            
            }   
            double simhitMaxE = 0.0;
            RawTrackerHit bestRawhit = null;
            for(RawTrackerHit rawhit : rawhitsonLayer){
                SimTrackerHit simhit = this.getBestSimHit(rawhit, rawtomc);
                if (simhit != null && simhit.getMCParticle() != null) {
                    double simhitEnergy = simhit.getdEdx();
                    if(simhitEnergy > simhitMaxE){
                        simhitMaxE = simhitEnergy;
                        bestRawhit = rawhit;
                    }
                }
            }
            rawhitsPerLayer.add(bestRawhit);
        }
        return rawhitsPerLayer;

    }   
        

    public MCParticle newgetTrackMCP(EventHeader event, Track track, String trackCollectionName, Boolean hitRequirement){

        Map <MCParticle, int[]>mcParticleMultiplicity = new HashMap<MCParticle, int[]>();

        //Retrieve rawhits to mc
        RelationalTable rawtomc = new BaseRelationalTable(RelationalTable.Mode.MANY_TO_MANY, RelationalTable.Weighting.UNWEIGHTED);
        if (event.hasCollection(LCRelation.class, "SVTTrueHitRelations")) {
            List<LCRelation> trueHitRelations = event.get(LCRelation.class, "SVTTrueHitRelations");
            for (LCRelation relation : trueHitRelations)
                if (relation != null && relation.getFrom() != null && relation.getTo() != null)
                    rawtomc.add(relation.getFrom(), relation.getTo());
        }

        MCParticle particle;
        for(TrackerHit hit : track.getTrackerHits()){
            Set<RawTrackerHit> rawhits = getRawHitsFromTrackerHit(hit, rawtomc);
            for(RawTrackerHit rawhit : rawhits){
                SimTrackerHit simhit = getBestSimHit(rawhit, rawtomc);
                if(simhit != null && simhit.getMCParticle() != null){
                    particle = simhit.getMCParticle();
                    if(!mcParticleMultiplicity.containsKey(particle)){
                        mcParticleMultiplicity.put(particle, new int[1]);
                        mcParticleMultiplicity.get(particle)[0] = 0;
                    }

                    mcParticleMultiplicity.get(particle)[0]++;
                }
            }

        }

        // Look for the MC particle that occurs the most of the track
        int maxValue = 0;
        particle = null;
        for(Map.Entry<MCParticle, int[]> entry : mcParticleMultiplicity.entrySet()){
            if(maxValue < entry.getValue()[0]){
                particle = entry.getKey();
                maxValue = entry.getValue()[0];
            }
        }

        int charge = -1* (int)Math.signum(track.getTrackStates().get(0).getOmega());
        if(charge > 0 && particle.getPDGID() == 11)
            return null;

        if(charge < 0 && particle.getPDGID() == -11)
            return null;

        if(charge < 0){
            plots1D.get(String.format("%s_ele_track_maxMCPmultiplicity",this.trackCollectionName)).fill(maxValue);
        }
        else{
            plots1D.get(String.format("%s_pos_track_maxMCPmultiplicity",this.trackCollectionName)).fill(maxValue);
        }

        if(hitRequirement == true){

            if(trackCollectionName.contains("GBLTracks")){
                if(maxValue > 9)
                    return particle;
                else
                    return null;
            }
            else if(trackCollectionName.contains("KalmanFullTracks")){
                if(maxValue > 9)
                    return particle;
                else
                    return null;
            }
            else
                return null;
        }
        else
            return particle;
    }

    private MCParticle getTrackMCP(Track track, EventHeader event){
        
        Map <MCParticle, int[]>mcParticleMultiplicity = new HashMap<MCParticle, int[]>();

        //Retrieve rawhits to mc
        RelationalTable rawtomc = new BaseRelationalTable(RelationalTable.Mode.MANY_TO_MANY, RelationalTable.Weighting.UNWEIGHTED);
        if (event.hasCollection(LCRelation.class, "SVTTrueHitRelations")) {
            List<LCRelation> trueHitRelations = event.get(LCRelation.class, "SVTTrueHitRelations");
            for (LCRelation relation : trueHitRelations)
                if (relation != null && relation.getFrom() != null && relation.getTo() != null)
                    rawtomc.add(relation.getFrom(), relation.getTo());
        }

        //Retrieve all simulated hits
        String MCHitInputCollectionName = "TrackerHits";
        List<SimTrackerHit> allsimhits = event.get(SimTrackerHit.class, MCHitInputCollectionName);

        MCParticle particle;
        for(TrackerHit hit : track.getTrackerHits()){
            List<RawTrackerHit> rawhits = hit.getRawHits();
            Set<RawTrackerHit> EmaxRawHitPerLayer = new HashSet<RawTrackerHit>();
            Map<RawTrackerHit,SimTrackerHit> rawToSimMap = new HashMap<RawTrackerHit, SimTrackerHit>();

            for(RawTrackerHit rawhit : rawhits){
                if(rawhit == null)
                    continue;
                Set<SimTrackerHit> simhits = rawtomc.allFrom(rawhit);
                SimTrackerHit rawhitSimhit = null;
                double simhitMaxE = 0.0;
                for(SimTrackerHit simhit : simhits){
                    if (simhit != null && simhit.getMCParticle() != null) {
                        double simhitEnergy = simhit.getdEdx();
                        if(simhitEnergy > simhitMaxE){
                            simhitMaxE = simhitEnergy;
                            rawhitSimhit = simhit;
                        }
                    }
                }
                if(rawhitSimhit != null)
                    rawToSimMap.put(rawhit, rawhitSimhit);
            }

            Map<RawTrackerHit, Integer> layerMap = new HashMap<RawTrackerHit, Integer>(); 
            Set<Integer> layers = new HashSet<Integer>();
            for(RawTrackerHit rawhit : rawToSimMap.keySet()){
                int layer = rawToSimMap.get(rawhit).getLayer(); 
                layerMap.put(rawhit, layer);
                layers.add(layer);
            }

            for(int layer : layers){
                List<RawTrackerHit> rawhitsonLayer = new ArrayList<RawTrackerHit>();
                for(Map.Entry<RawTrackerHit,Integer> entry : layerMap.entrySet()){
                    if(entry.getValue().equals(layer)){
                        rawhitsonLayer.add(entry.getKey());
                    }
                
                }   
                double simhitMaxE = 0.0;
                SimTrackerHit bestSimhit = null;
                for(RawTrackerHit rawhit : rawhitsonLayer){
                    SimTrackerHit simhit = rawToSimMap.get(rawhit);
                    if (simhit != null && simhit.getMCParticle() != null) {
                        double simhitEnergy = simhit.getdEdx();
                        if(simhitEnergy > simhitMaxE){
                            simhitMaxE = simhitEnergy;
                            bestSimhit = simhit;
                        }
                    }
                }

                if(bestSimhit != null && bestSimhit.getMCParticle() != null){
                    particle = bestSimhit.getMCParticle();
                    if(!mcParticleMultiplicity.containsKey(particle)){
                        mcParticleMultiplicity.put(particle, new int[1]);
                        mcParticleMultiplicity.get(particle)[0] = 0;
                    }

                    mcParticleMultiplicity.get(particle)[0]++;
                }
            }
            



        }

        // Look for the MC particle that occurs the most of the track
        int maxValue = 0;
        int minhits = 9;
        particle = null;
        for(Map.Entry<MCParticle, int[]> entry : mcParticleMultiplicity.entrySet()){
            if(maxValue < entry.getValue()[0]){
                particle = entry.getKey();
                maxValue = entry.getValue()[0];
            }
        }

        int charge = -1* (int)Math.signum(track.getTrackStates().get(0).getOmega());
        if(charge < 0){
            plots1D.get(String.format("%s_ele_track_maxMCPmultiplicity",this.trackCollectionName)).fill(maxValue);
        }
        else{
            plots1D.get(String.format("%s_pos_track_maxMCPmultiplicity",this.trackCollectionName)).fill(maxValue);
        }
        if(this.trackCollectionName.contains("GBLTracks")){
            if(maxValue > minhits)
                return particle;
            else
                return null;
        }
        else if(this.trackCollectionName.contains("KalmanFullTracks")){
            if(maxValue > minhits)
                return particle;
            else
                return null;
        }
        else
            return null;
    }

}


