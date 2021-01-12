package org.hps.recon.tracking.kalman;

import hep.aida.IAnalysisFactory;
import hep.aida.IHistogram1D;
import hep.aida.IHistogramFactory;
import hep.aida.ITree;
import hep.aida.ref.rootwriter.RootFileStore;

import org.hps.recon.ecal.cluster.ClusterUtilities;
import org.hps.recon.tracking.TrackUtils;
import org.hps.recon.tracking.TrackData;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.List;
import java.util.ArrayList;

import org.lcsim.event.Track;
import org.lcsim.event.TrackState;
import org.lcsim.event.Cluster;
import org.lcsim.event.RelationalTable;
import org.lcsim.event.EventHeader;
import org.lcsim.event.base.BaseRelationalTable;
import org.lcsim.event.LCRelation;

public class TrackClusterMatcher2019 {


    //Used for plots
    boolean enablePlots = false;
    private ITree tree;
    private IHistogramFactory histogramFactory;
    private Map<String, IHistogram1D> plots1D;
    private String trackType;

    RelationalTable hitToRotated = null;
    RelationalTable hitToStrips = null;

    public TrackClusterMatcher2019(String trackCollectionName) {
        trackType = trackCollectionName;
        System.out.println("Matching " + trackType + "to Ecal Clusters");
    }

    public void enablePlots(boolean enablePlots) {
        this.enablePlots = enablePlots;
        if (enablePlots ==true) {
            this.bookHistograms();
        }
    }

    public void saveHistograms() {
        System.out.println("Saving Histogram for " + this.trackType);
        String rootFile = String.format("%s_TrackClusterMatching.root",this.trackType);
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

        System.out.println("BOOKING HISTOGRAMS FOR " +  this.trackType);
        String trackType = this.trackType;
        plots1D = new HashMap<String, IHistogram1D>();
        tree = IAnalysisFactory.create().createTreeFactory().create();
        histogramFactory = IAnalysisFactory.create().createHistogramFactory(tree);

        //Timing Plots
        plots1D.put(String.format("%s_ElectronTrack-Cluster_dt",trackType), histogramFactory.createHistogram1D(String.format("%s_ElectronTrack-Cluster_dt",trackType),  200, -200, 200));
        plots1D.put(String.format("%s_PositronTrack-Cluster_dt",trackType), histogramFactory.createHistogram1D(String.format("%s_PositronTrack-Cluster_dt",trackType),  200, -200, 200));
        plots1D.put(String.format("%s_PositronTrackTime",trackType), histogramFactory.createHistogram1D(String.format("%s_PositronTrackTime",trackType),  200, -200, 200));
        plots1D.put(String.format("%s_ElectronTrackTime",trackType), histogramFactory.createHistogram1D(String.format("%s_ElectronTrackTime",trackType),  200, -200, 200));
        plots1D.put(String.format("%s_Cluster_Timing_(woffset)",trackType), histogramFactory.createHistogram1D(String.format("%s_Cluster_Timing_(woffset)",trackType),  200, -200, 200));

        //Track to Ecal Extrapolation
        plots1D.put(String.format("%s_PositronTrack@ECal_ECalCluster_dx",trackType),histogramFactory.createHistogram1D(String.format("%s_PositronTrack@ECal_ECalCluster_dx",trackType),400, -200, 200));
        plots1D.put(String.format("%s_PositronTrack@ECal_ECalCluster_dy",trackType),histogramFactory.createHistogram1D(String.format( "%s_PositronTrack@ECal_ECalCluster_dy",trackType),400, -200, 200));
        plots1D.put(String.format("%s_PositronTrack@ECal_ECalCluster_dz",trackType),histogramFactory.createHistogram1D(String.format( "%s_PositronTrack@ECal_ECalCluster_dz",trackType),100, -50,50));
        plots1D.put(String.format("%s_PositronTrack@ECal_ECalCluster_dr",trackType),histogramFactory.createHistogram1D(String.format( "%s_PositronTrack@ECal_ECalCluster_dr",trackType),100, -50, 150));
        plots1D.put(String.format("%s_ElectronTrack@ECal_ECalCluster_dx",trackType),histogramFactory.createHistogram1D(String.format("%s_ElectronTrack@ECal_ECalCluster_dx",trackType),400, -200, 200));
        plots1D.put(String.format("%s_ElectronTrack@ECal_ECalCluster_dy",trackType),histogramFactory.createHistogram1D(String.format( "%s_ElectronTrack@ECal_ECalCluster_dy",trackType),400, -200, 200));
        plots1D.put(String.format("%s_ElectronTrack@ECal_ECalCluster_dz",trackType),histogramFactory.createHistogram1D(String.format( "%s_ElectronTrack@ECal_ECalCluster_dz",trackType),100, -50,50));
        plots1D.put(String.format("%s_ElectronTrack@ECal_ECalCluster_dr",trackType),histogramFactory.createHistogram1D(String.format( "%s_ElectronTrack@ECal_ECalCluster_dr",trackType),100, -50,150));

        // Energy/Momentum plots
        plots1D.put(String.format("%s_ele_Track_Cluster_EdivP",trackType), histogramFactory.createHistogram1D(String.format("%s_ele_Track_Cluster_EdivP",trackType),  200, -10, 10));
        plots1D.put(String.format("%s_pos_Track_Cluster_EdivP",trackType), histogramFactory.createHistogram1D(String.format("%s_pos_Track_Cluster_EdivP",trackType),  200, -10, 10));
        plots1D.put(String.format("%s_ele_Track_Momentum",trackType), histogramFactory.createHistogram1D(String.format("%s_ele_Track_Momentum",trackType),  100, 0, 5));
        plots1D.put(String.format("%s_pos_Track_Momentum",trackType), histogramFactory.createHistogram1D(String.format("%s_pos_Track_Momentum",trackType),  100, 0, 5));
        plots1D.put(String.format("%s_cluster_energy",trackType), histogramFactory.createHistogram1D(String.format("%s_cluster_energy",trackType),  100, 0, 5));
    }

    public void plotClusterEnergy(EventHeader event, String tag, String ecalClustersCollectionName){
        System.out.println("[TrackClusterMatcher2019] Event Number: " + event.getEventNumber());

        List<Cluster> clusters = event.get(Cluster.class, ecalClustersCollectionName);    
        for(Cluster cluster : clusters) {
            double clusterEnergy = cluster.getEnergy();
            plots1D.get(String.format("%s_cluster_energy_fromEvent_%s",trackType,tag)).fill(clusterEnergy);
            System.out.println("[TrackClusterMatcher2019] " + this.trackType + " Cluster Energy: " + clusterEnergy);
        }
    }

    public void plotClusterEnergy(List<Cluster> clusters, String tag){
        for(Cluster cluster : clusters) {
            double clusterEnergy = cluster.getEnergy();
            plots1D.get(String.format("%s_cluster_energy_fromList_%s",trackType,tag)).fill(clusterEnergy);
            System.out.println("List Cluster Energy: " + clusterEnergy);
        }
    }

    public Map<Track,Cluster> trackClusterMatcher(List<Track> tracks, EventHeader event,  String trackCollectionName, List<Cluster> clusters, double trackClusterTimeOffset)  {

        //Input collection of Tracks, with trackCollectionName, and collection
        //of Ecal Clusters
        //Method matches unique Ecal clusters to Tracks based on closest
        //distance, within a specific time window
        //Output is a map between Tracks and matched Cluster
        //If no cluster is matched to a Track, Map contains Track + Null
        //cluster

        
        //Map of position residuals between all track+cluster combinations
        if(tracks == null || tracks.isEmpty()){
            System.out.println("Track list given to KFTrackEcalClusterMatcher is Empty!");
            return null;
        }

        Map<Track, Map<Cluster, Double>> trackClusterResMap = new HashMap<Track, Map<Cluster, Double>>(); 

        //First, gather all necessary Track information
        String trackType = trackCollectionName;

        //Relation Table required to retrieve kalman track time through
        //TrackData class
        hitToRotated = TrackUtils.getHitToRotatedTable(event);
        hitToStrips = TrackUtils.getHitToStripsTable(event);
        List<TrackData> TrackData;
        RelationalTable trackToData = new BaseRelationalTable(RelationalTable.Mode.ONE_TO_ONE, RelationalTable.Weighting.UNWEIGHTED);
        List<LCRelation> trackRelations;
        TrackData trackdata;
        if (trackCollectionName.contains("KalmanFullTracks")) {
            TrackData = event.get(TrackData.class, "KFTrackData");
            trackRelations = event.get(LCRelation.class, "KFTrackDataRelations");
            for (LCRelation relation : trackRelations) {
                if (relation != null && relation.getTo() != null){
                    trackToData.add(relation.getFrom(), relation.getTo());
                }
            }
        }


        for(Track track : tracks) {

            //charge sign must be flipped by factor of -1 (WHY!?)
            int charge = -1* (int) Math.signum(track.getTrackStates().get(0).getOmega());
            double trackt;
            double tracktOffset = 4; //Track time distribution is centered on -4ns, added offset to center ~0

            double trackx;
            double tracky;
            double trackz;
            double dxoffset; //Track x-position at Ecal distriution not centered on 0. Offset varies depending on charge

            //Track parameters that are useful for implementing cuts on
            //potential cluster matches
            double tanlambda = track.getTrackParameter(4);

            if (trackType.contains("GBLTracks")){
                trackt = TrackUtils.getTrackTime(track, hitToStrips, hitToRotated);
                trackx = TrackUtils.getTrackStateAtECal(track).getReferencePoint()[1]; 
                tracky = TrackUtils.getTrackStateAtECal(track).getReferencePoint()[2];
                trackz = TrackUtils.getTrackStateAtECal(track).getReferencePoint()[0];
                //trackP = track.getTrackStates().get(0).getMomentum(); 

                //electron GBLTracks show x position bias at +5.5 mm. Offset
                //accounted for below (IS THIS OKAY TO HARDCODE?)
                if(charge < 0)
                    dxoffset = -5.5;
                else
                    dxoffset = 0.0;
            }

            //KFTracks
            else {
                trackdata = (TrackData) trackToData.from(track);
                trackt = trackdata.getTrackTime();
                //KF TrackState at ecal stored as the last TrackState in
                //KalmanInterface.java
                TrackState ts_ecal = track.getTrackStates().get(track.getTrackStates().size()-1); //Be careful about the coordinate frame used for this track state. It is different between current master and pass1-dev-fix branches.
                //trackP = track.getTrackStates().get(0).getMomentum(); 
                double[] ts_ecalPos = ts_ecal.getReferencePoint();
                if(ts_ecalPos[0] > 1000){
                    trackx = ts_ecalPos[1];
                    tracky = ts_ecalPos[2];
                    trackz = ts_ecalPos[0];
                }
                else{
                    trackx = ts_ecalPos[0];
                    tracky = ts_ecalPos[1];
                    trackz = ts_ecalPos[2];
                }
                if(charge < 0)
                    dxoffset = 3.3; //KF ele tracks have x-position bias of -3.3 mm, hardcode offset + 3.3
                else
                    dxoffset = - 3.6; //Similar case as above 
            }

            //Track momentum magnitude
            double[] trackP = track.getMomentum();
            double trackPmag = Math.sqrt(Math.pow(trackP[0],2) + Math.pow(trackP[1],2) + Math.pow(trackP[2],2));
            //FEE's Only Cut
            if(trackPmag < 2.0 || trackPmag > 2.5)
                continue;

            //Plots
            if(enablePlots){
                if (charge > 0) {
                    plots1D.get(String.format("%s_PositronTrackTime",trackType)).fill(trackt);
                    plots1D.get(String.format("%s_pos_Track_Momentum",trackType)).fill(trackPmag);
                }
                else {
                    plots1D.get(String.format("%s_ElectronTrackTime",trackType)).fill(trackt);
                    plots1D.get(String.format("%s_ele_Track_Momentum",trackType)).fill(trackPmag);
                }
            }

            //Begin Cluster Matching Algorithm
            Map<Cluster, Double> clusterResMap = new HashMap<Cluster, Double>();
            Cluster matchedCluster = null;

            double smallestdt = Double.MAX_VALUE;
            double smallestdr = Double.MAX_VALUE;
            //define time and position cuts for Track-Cluster matching,
            //determined by distributions
            double tcut = 4.0;
            double xcut = 10.0;
            double ycut = 10.0;

            //Loop over all clusters, looking for best match to current track
            for(Cluster cluster : clusters) {
                double clusterEnergy = cluster.getEnergy();
                double clusTime = ClusterUtilities.getSeedHitTime(cluster);
                double dt = clusTime - trackClusterTimeOffset - trackt + tracktOffset;

                double clusterx = cluster.getPosition()[0];
                double clustery = cluster.getPosition()[1];
                double clusterz = cluster.getPosition()[2];
                double dx = clusterx - trackx + dxoffset;
                double dy = clustery - tracky;
                double dz = clusterz - trackz;
                double dr = Math.sqrt(Math.pow(clusterx-trackx,2) + Math.pow(clustery-tracky,2));

                /*
                //Ecal fiducial cuts
                if(clusterx < 0 && charge > 0)
                    continue;
                if(clusterx > 0 && charge < 0)
                    continue;
                if(clustery > 0 && tanlambda < 0)
                    continue;
                if(clustery < 0 && tanlambda > 0)
                    continue;
                    */

                //Plot of cluster energy / track momentum
                if(enablePlots){
                    plots1D.get(String.format("%s_cluster_energy",trackType)).fill(clusterEnergy);
                    if(charge < 0)
                        plots1D.get(String.format("%s_ele_Track_Cluster_EdivP",trackType)).fill(clusterEnergy/trackPmag);
                    else
                        plots1D.get(String.format("%s_pos_Track_Cluster_EdivP",trackType)).fill(clusterEnergy/trackPmag);
                }

                //If position and time residual cuts are passed, build map of
                //all cluster position residuals with this track
                if((Math.abs(dt) < tcut) && (Math.abs(dx) < xcut) && (Math.abs(dy) < ycut) ) {
                    clusterResMap.put(cluster, dr);
                }

                if(enablePlots) {
                    plots1D.get(String.format("%s_Cluster_Timing_(woffset)",trackType)).fill(clusTime - trackClusterTimeOffset);
                    if((Math.abs(dt) < tcut) && (Math.abs(dx) < xcut) && (Math.abs(dy) < ycut) ) {
                        if(charge > 0) {
                            //Time residual plot
                            plots1D.get(String.format("%s_PositronTrack-Cluster_dt",trackType)).fill(dt);

                            //Kalman Extrapolated Residuals
                            plots1D.get(String.format("%s_PositronTrack@ECal_ECalCluster_dx",trackType)).fill(dx);
                            plots1D.get(String.format("%s_PositronTrack@ECal_ECalCluster_dy",trackType)).fill(dy);
                            plots1D.get(String.format("%s_PositronTrack@ECal_ECalCluster_dz",trackType)).fill(dz);
                            plots1D.get(String.format("%s_PositronTrack@ECal_ECalCluster_dr",trackType)).fill(dr);
                        }
                        else {

                            //Time residual plot
                            plots1D.get(String.format("%s_ElectronTrack-Cluster_dt",trackType)).fill(dt);

                            //Kalman Extrapolated Residuals
                            plots1D.get(String.format("%s_ElectronTrack@ECal_ECalCluster_dx",trackType)).fill(dx);
                            plots1D.get(String.format("%s_ElectronTrack@ECal_ECalCluster_dy",trackType)).fill(dy);
                            plots1D.get(String.format("%s_ElectronTrack@ECal_ECalCluster_dz",trackType)).fill(dz);
                            plots1D.get(String.format("%s_ElectronTrack@ECal_ECalCluster_dr",trackType)).fill(dr);
                        }
                    }
                }

            }

            //trackClusterResMap is a map of Tracks to the position residual of
            //each potential cluster match
            trackClusterResMap.put(track, clusterResMap);
        }

        //Given the mapping between all Tracks, and all potential cluster
        //matches, match Tracks to the Clusters that are closest in position
        //Algorithm checks for clusters matched to multiple Tracks, and sorts
        //them until only unique matches exist

        //trackMinResClusterMap maps Track to the Cluster with the smallest dr
        Map<Track,Cluster> trackMinResClusterMap = new HashMap<Track, Cluster>();

        //check map for the same Cluster being matched to multiple Tracks
        //If found, keep Track-Cluster match with smallest dr
        //Repeat matching process
        for(int i=0; i < clusters.size(); i++){
            trackMinResClusterMap = getTrackMinResClusterMap(trackClusterResMap);
            trackClusterResMap = checkDuplicateClusterMatching(trackClusterResMap,trackMinResClusterMap);
        }
        trackMinResClusterMap = getTrackMinResClusterMap(trackClusterResMap);
        return trackMinResClusterMap;
    }


    public Map<Track, Map<Cluster,Double>> checkDuplicateClusterMatching(Map<Track, Map<Cluster,Double>> trackClusterResMap, Map<Track, Cluster> trackMinResClusterMap){
        
        boolean duplicateCluster = false;
        List<Track> sharedClusterTracks = new ArrayList<Track>();
        List<Track> skipTracks = new ArrayList<Track>();
        
        for(Track track : trackMinResClusterMap.keySet()){
            Map<Track, Cluster> trackswDuplicateClusters = new HashMap<Track, Cluster>();
            if(skipTracks.contains(track))
                continue;
            Cluster smallestdrCluster = trackMinResClusterMap.get(track);
            if(smallestdrCluster == null)
                continue;
            for(Track otherTrack : trackMinResClusterMap.keySet()){
                if(skipTracks.contains(track))
                    continue;
                if(otherTrack == track)
                    continue;
                Cluster othersmallestdrCluster = trackMinResClusterMap.get(otherTrack);
                if(othersmallestdrCluster == smallestdrCluster)
                {
                    duplicateCluster = true;
                    trackswDuplicateClusters.put(track, smallestdrCluster);
                    trackswDuplicateClusters.put(otherTrack, othersmallestdrCluster);
                }
            }

            double smallestdr = 99999.0;
            Track smallestdrTrack = null;
            if(trackswDuplicateClusters == null)
                return trackClusterResMap;
            for(Track duptrack : trackswDuplicateClusters.keySet()){
                double dr = trackClusterResMap.get(duptrack).get(trackswDuplicateClusters.get(duptrack));
                if(dr < smallestdr){
                    smallestdr = dr;
                    smallestdrTrack = duptrack;
                }
            }
            for(Track duptrack : trackswDuplicateClusters.keySet()){
                skipTracks.add(duptrack);
                if(duptrack != smallestdrTrack){
                    trackClusterResMap.get(duptrack).remove(trackswDuplicateClusters.get(duptrack));
                }
            }
        }
        return trackClusterResMap;
    }

    public Map<Track,Cluster>  getTrackMinResClusterMap(Map<Track, Map<Cluster, Double>> trackClusterResMap){

        //inputs a mapping of tracks with residuals for each possible cluster
        //from all clusters in the map, for each track, match the cluster with
        //the smallest position residual to that track
        //build output map of track -> closest cluster
        Map<Track,Cluster> Map = new HashMap<Track, Cluster>();
        for(Track track : trackClusterResMap.keySet()){
            double smallestdr = 99999.0;
            Cluster smallestdrCluster = null;
            Map<Cluster, Double> clusterResMap = trackClusterResMap.get(track);
            for(Cluster c : clusterResMap.keySet()){
                double dr = clusterResMap.get(c);
                if(dr < smallestdr){
                    smallestdr = dr;
                    smallestdrCluster = c;
                }
            }
            Map.put(track, smallestdrCluster);
        }
        return Map;
        
    }


}
