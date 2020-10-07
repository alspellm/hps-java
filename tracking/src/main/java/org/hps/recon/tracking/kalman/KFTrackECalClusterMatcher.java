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

public class KFTrackECalClusterMatcher {


    //Used for plots
    boolean enablePlots;
    private ITree tree;
    private IHistogramFactory histogramFactory;
    private Map<String, IHistogram1D> plots1D;
    private String trackType;

    RelationalTable hitToRotated = null;
    RelationalTable hitToStrips = null;

    public KFTrackECalClusterMatcher(String trackCollectionName) {
        trackType = trackCollectionName;
        System.out.println("TRACK TYPE IS " + trackType);
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

       /* 
        //RK Extrapolated Residuals
        plots1D.put(String.format("RK_PositronTrack@ECal_ECalCluster_dx",histogramFactory.createHistogram1D(String.format("RK_PositronTrack@ECal_ECalCluster_dx",400, -200, 200));
        plots1D.put(String.format("RK_PositronTrack@ECal_ECalCluster_dy",histogramFactory.createHistogram1D(String.format( "RK_PositronTrack@ECal_ECalCluster_dy",400, -200, 200));
        plots1D.put(String.format("RK_PositronTrack@ECal_ECalCluster_dz",histogramFactory.createHistogram1D(String.format( "RK_PositronTrack@ECal_ECalCluster_dz",100, -50,50));
        plots1D.put(String.format("RK_ElectronTrack@ECal_ECalCluster_dx",histogramFactory.createHistogram1D(String.format("RK_ElectronTrack@ECal_ECalCluster_dx",400, -200, 200));
        plots1D.put(String.format("RK_ElectronTrack@ECal_ECalCluster_dy",histogramFactory.createHistogram1D(String.format( "RK_ElectronTrack@ECal_ECalCluster_dy",400, -200, 200));
        plots1D.put(String.format("RK_ElectronTrack@ECal_ECalCluster_dz",histogramFactory.createHistogram1D(String.format( "RK_ElectronTrack@ECal_ECalCluster_dz",100, -50,50));
    */
    }
    public Cluster trackClusterMatcher(Track TrackHPS, String trackCollectionName, int Charge, List<Cluster> Clusters, double trackTime, double trackClusterTimeOffset)  {

        //HPSTrack
        Track track = TrackHPS;
        int charge = Charge;
        double tanlambda = track.getTrackParameter(4);
        String trackType = trackCollectionName;
        double tracktOffset = 4; //track time distributions show mean at -4 ns
        double trackt = trackTime;
        double trackx;
        double tracky;
        double trackz;
        double dxoffset;
        List<Cluster> clusters = Clusters;

        if(trackType.contains("GBLTracks")) {
            trackx = TrackUtils.getTrackStateAtECal(track).getReferencePoint()[1]; 
            tracky = TrackUtils.getTrackStateAtECal(track).getReferencePoint()[2];
            trackz = TrackUtils.getTrackStateAtECal(track).getReferencePoint()[0];
            if(charge < 0)
                dxoffset = -5.5;
            else
                dxoffset = 0.0;
        }

        else {
            TrackState ts_ecal = track.getTrackStates().get(track.getTrackStates().size()-1);
            double[] ts_ecalPos = ts_ecal.getReferencePoint();
            trackx = ts_ecalPos[0];
            tracky = ts_ecalPos[1];
            trackz = ts_ecalPos[2];
            dxoffset = 0.0;
        }


        /*
        //Track state at ecal via RK extrap
        TrackState ts_ecalRK = track.getTrackStates().get(track.getTrackStates().size()-2);
        Hep3Vector ts_ecalPos_RK = new BasicHep3Vector(ts_ecalRK.getReferencePoint());
        ts_ecalPos_RK = CoordinateTransformations.transformVectorToDetector(ts_ecalPos_RK);
        */


        Cluster matchedCluster = null;
        double smallestdt = Double.MAX_VALUE;
        double smallestdr = Double.MAX_VALUE;
        double offset = trackClusterTimeOffset;
        double tcut = 4.0;
        double xcut = 20.0;
        double ycut = 15.0;

        for(Cluster cluster : clusters) {
            double clusTime = ClusterUtilities.getSeedHitTime(cluster);
            double dt = clusTime - trackClusterTimeOffset - trackt + tracktOffset;

            double clusxpos = cluster.getPosition()[0];
            double clusypos = cluster.getPosition()[1];
            double cluszpos = cluster.getPosition()[2];
            double dx = clusxpos - trackx + dxoffset;
            double dy = clusypos - tracky;
            double dz = cluszpos - trackz;
            double dr = Math.sqrt(Math.pow(clusxpos-trackx,2) + Math.pow(clusypos-tracky,2));
            if(clusxpos < 0 && charge > 0)
                continue;
            if(clusxpos > 0 && charge < 0)
                continue;
            if(clusypos > 0 && tanlambda < 0)
                continue;
            if(clusypos < 0 && tanlambda > 0)
                continue;

            /*
            //via RK Extrap
            double dxRK = clusPos[0]-ts_ecalPos_RK.x();
            double dyRK = clusPos[1]-ts_ecalPos_RK.y();
            double dzRK = clusPos[2]-ts_ecalPos_RK.z();
            */


            //Extremely simplified track cluster matching. Cluster that passes
            //position cuts and has closest time is matched. This needs to be
            //updated to a real algorithm.
            if((Math.abs(dt) < tcut) && (Math.abs(dx) < xcut) && (Math.abs(dy) < ycut) ) {
                System.out.println("KF cluster passing selection found");
                if(Math.abs(dr) < smallestdr) {
                    smallestdr = Math.abs(dr);
                    matchedCluster = cluster;
                }
            }

            
            if(enablePlots) {
                System.out.println("Filling Histograms for " + trackType);
                plots1D.get(String.format("%s_Cluster_Timing_(woffset)",trackType)).fill(clusTime-offset);
                if((Math.abs(dt) < tcut) && (Math.abs(dx) < xcut) && (Math.abs(dy) < ycut) ) {
                    if(charge > 0) {
                        //Time residual plot
                        plots1D.get(String.format("%s_PositronTrack-Cluster_dt",trackType)).fill(dt);

                        //Kalman Extrapolated Residuals
                        plots1D.get(String.format("%s_PositronTrack@ECal_ECalCluster_dx",trackType)).fill(dx);
                        plots1D.get(String.format("%s_PositronTrack@ECal_ECalCluster_dy",trackType)).fill(dy);
                        plots1D.get(String.format("%s_PositronTrack@ECal_ECalCluster_dz",trackType)).fill(dz);
                        plots1D.get(String.format("%s_PositronTrack@ECal_ECalCluster_dr",trackType)).fill(dr);

                        /*
                        //RK Extrapolated Residuals
                        plots1D.get(String.format("RK_PositronTrack@ECal_ECalCluster_dx").fill(dxRK);
                        plots1D.get(String.format("RK_PositronTrack@ECal_ECalCluster_dy").fill(dyRK);
                        plots1D.get(String.format("RK_PositronTrack@ECal_ECalCluster_dz").fill(dzRK);
                        */
                    }
                    else {

                        //Time residual plot
                        plots1D.get(String.format("%s_ElectronTrack-Cluster_dt",trackType)).fill(dt);

                        //Kalman Extrapolated Residuals
                        plots1D.get(String.format("%s_ElectronTrack@ECal_ECalCluster_dx",trackType)).fill(dx);
                        plots1D.get(String.format("%s_ElectronTrack@ECal_ECalCluster_dy",trackType)).fill(dy);
                        plots1D.get(String.format("%s_ElectronTrack@ECal_ECalCluster_dz",trackType)).fill(dz);
                        plots1D.get(String.format("%s_ElectronTrack@ECal_ECalCluster_dr",trackType)).fill(dr);
                        /*
                        //RK Extrapolated Residuals
                        plots1D.get(String.format("RK_ElectronTrack@ECal_ECalCluster_dx").fill(dxRK);
                        plots1D.get(String.format("RK_ElectronTrack@ECal_ECalCluster_dy").fill(dyRK);
                        plots1D.get(String.format("RK_ElectronTrack@ECal_ECalCluster_dz").fill(dzRK);
                        */
                    }
                }
            }
            
        }

        if(enablePlots){
            if (charge > 0) {
                plots1D.get(String.format("%s_PositronTrackTime",trackType)).fill(trackt);
            }
            else {
                plots1D.get(String.format("%s_ElectronTrackTime",trackType)).fill(trackt);
            }
        }
        if(matchedCluster == null){
            System.out.println("No matching cluster found for KF track at ECal");
            return null;
        }
        else {
            return matchedCluster;
        }
    }

    public Map<Track,Cluster> newtrackClusterMatcher(List<Track> tracks, RelationalTable trackToData, RelationalTable hitToRotated, RelationalTable hitToStrips,  String trackCollectionName, List<Cluster> Clusters, double TrackClusterTimeOffset)  {

        String trackType = trackCollectionName;
        TrackData trackdata;
        Map<Track, Map<Cluster, Double>> trackClusterResMap = new HashMap<Track, Map<Cluster, Double>>(); 
        for(Track track : tracks) {
            int charge = -1* (int) Math.signum(track.getTrackStates().get(0).getOmega());
            double trackt;
            double tracktOffset = 4; //track time distributions show mean at -4 ns
            if (trackType.contains("GBLTracks")){
                trackt = TrackUtils.getTrackTime(track, hitToStrips, hitToRotated);
            }
            else {
                trackdata = (TrackData) trackToData.from(track);
                trackt = trackdata.getTrackTime();
            }

            double tanlambda = track.getTrackParameter(4);
            double trackx;
            double tracky;
            double trackz;
            double dxoffset;
            List<Cluster> clusters = Clusters;

            if(trackType.contains("GBLTracks")) {
                trackx = TrackUtils.getTrackStateAtECal(track).getReferencePoint()[1]; 
                tracky = TrackUtils.getTrackStateAtECal(track).getReferencePoint()[2];
                trackz = TrackUtils.getTrackStateAtECal(track).getReferencePoint()[0];
                if(charge < 0)
                    dxoffset = -5.5;
                else
                    dxoffset = 0.0;
            }

            else {
                TrackState ts_ecal = track.getTrackStates().get(track.getTrackStates().size()-1);
                double[] ts_ecalPos = ts_ecal.getReferencePoint();
                trackx = ts_ecalPos[0];
                tracky = ts_ecalPos[1];
                trackz = ts_ecalPos[2];
                dxoffset = 0.0;
            }


            /*
            //Track state at ecal via RK extrap
            TrackState ts_ecalRK = track.getTrackStates().get(track.getTrackStates().size()-2);
            Hep3Vector ts_ecalPos_RK = new BasicHep3Vector(ts_ecalRK.getReferencePoint());
            ts_ecalPos_RK = CoordinateTransformations.transformVectorToDetector(ts_ecalPos_RK);
            */

            Cluster matchedCluster = null;
            double smallestdt = Double.MAX_VALUE;
            double smallestdr = Double.MAX_VALUE;
            double trackClusterTimeOffset = TrackClusterTimeOffset;
            //double tcut = 4.0;
            //double xcut = 20.0;
            //double ycut = 15.0;
            double tcut = 999;
            double xcut = 999;
            double ycut = 999;

            Map<Cluster, Double> clusterResMap = new HashMap<Cluster, Double>();
            for(Cluster cluster : clusters) {
                double clusTime = ClusterUtilities.getSeedHitTime(cluster);
                double dt = clusTime - trackClusterTimeOffset - trackt + tracktOffset;

                double clusterx = cluster.getPosition()[0];
                double clustery = cluster.getPosition()[1];
                double clusterz = cluster.getPosition()[2];
                double dx = clusterx - trackx + dxoffset;
                double dy = clustery - tracky;
                double dz = clusterz - trackz;
                double dr = Math.sqrt(Math.pow(clusterx-trackx,2) + Math.pow(clustery-tracky,2));
                //if(clusterx < 0 && charge > 0)
                  //  continue;
                //if(clusterx > 0 && charge < 0)
                  //  continue;
                //if(clustery > 0 && tanlambda < 0)
                  //  continue;
                //if(clustery < 0 && tanlambda > 0)
                  //  continue;
                if((Math.abs(dt) < tcut) && (Math.abs(dx) < xcut) && (Math.abs(dy) < ycut) ) {
                    clusterResMap.put(cluster, dr);
                }



                /*
                //via RK Extrap
                double dxRK = clusPos[0]-ts_ecalPos_RK.x();
                double dyRK = clusPos[1]-ts_ecalPos_RK.y();
                double dzRK = clusPos[2]-ts_ecalPos_RK.z();
                */


                //Extremely simplified track cluster matching. Cluster that passes
                //position cuts and has closest time is matched. This needs to be
                //updated to a real algorithm.
                /*
                if((Math.abs(dt) < tcut) && (Math.abs(dx) < xcut) && (Math.abs(dy) < ycut) ) {
                    System.out.println("KF cluster passing selection found");
                    if(Math.abs(dr) < smallestdr) {
                        smallestdr = Math.abs(dr);
                        matchedCluster = cluster;
                    }
                }
                */
            }
            trackClusterResMap.put(track, clusterResMap);
            
        }
        System.out.println("Track length: " + tracks.size());
        System.out.println("Clusters length: " + Clusters.size());

        Map<Track,Cluster> trackMinResClusterMap = new HashMap<Track, Cluster>();
        for(int i=0; i < Clusters.size(); i++){
            trackMinResClusterMap = getTrackMinResClusterMap(trackClusterResMap);
            trackClusterResMap = checkDuplicateClusterMatching(trackClusterResMap,trackMinResClusterMap);
        }
        trackMinResClusterMap = getTrackMinResClusterMap(trackClusterResMap);
        return trackMinResClusterMap;

        //Map<Track,Cluster> matchedTrackClusterMap = new HashMap<Track, Cluster>();
        //for(Track track : trackMinResClusterMap.keySet()){
          // matchedTrackClusterMap.put(track,trackMinResClusterMap.get(track)); 
       // }


    }

    public void testList(List<Integer> list){
        list.remove(1);
        list.remove(3);

    }

    public Map<Track, Map<Cluster,Double>> checkDuplicateClusterMatching(Map<Track, Map<Cluster,Double>> trackClusterResMap, Map<Track, Cluster> trackMinResClusterMap){
        
        boolean duplicateCluster = false;
        List<Track> sharedClusterTracks = new ArrayList<Track>();
        List<Track> skipTracks = new ArrayList<Track>();
        
        for(Track track : trackMinResClusterMap.keySet()){
            Map<Track, Cluster> trackswDuplicateClusters = new HashMap<Track, Cluster>();
            if(skipTracks.contains(track))
                continue;
            System.out.println("Checking track for shared matching clusters");
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
                    System.out.println("Track matched to same Cluster dr: "+ trackClusterResMap.get(track).get(smallestdrCluster));
                    System.out.println("Track matched to same Cluster dr: "+ trackClusterResMap.get(otherTrack).get(othersmallestdrCluster));
                }
            }

            double smallestdr = 99999.0;
            Track smallestdrTrack = null;
            if(trackswDuplicateClusters == null)
                return trackClusterResMap;
            System.out.println("Size of tracksDuplicateCluster: " + trackswDuplicateClusters.size());
            for(Track duptrack : trackswDuplicateClusters.keySet()){
                double dr = trackClusterResMap.get(duptrack).get(trackswDuplicateClusters.get(duptrack));
                System.out.println("checking repeat list for dr value: " + dr);
                if(dr < smallestdr){
                    smallestdr = dr;
                    smallestdrTrack = duptrack;
                }
            }
            System.out.println("Duplicate Track with the smallest dr: " + smallestdr);
            for(Track duptrack : trackswDuplicateClusters.keySet()){
                skipTracks.add(duptrack);
                if(duptrack != smallestdrTrack){
                    trackClusterResMap.get(duptrack).remove(trackswDuplicateClusters.get(duptrack));
                }
            }
        }
        System.out.println("Finished checking for shared matching clusters");
        return trackClusterResMap;
    }

    public Map<Track,Cluster>  getTrackMinResClusterMap(Map<Track, Map<Cluster, Double>> trackClusterResMap){

        Map<Track,Cluster> trackMinResClusterMap = new HashMap<Track, Cluster>();
        for(Track track : trackClusterResMap.keySet()){
            System.out.println("Matching Track to Cluster with minimum dr");
            double smallestdr = 99999.0;
            Cluster smallestdrCluster = null;
            Map<Cluster, Double> clusterResMap = trackClusterResMap.get(track);
            for(Cluster c : clusterResMap.keySet()){
                double dr = clusterResMap.get(c);
                if(dr < smallestdr){
                    smallestdr = dr;
                    smallestdrCluster = c;
                }
                System.out.println("[KFTEM] track dr: " + dr);
            }
            trackMinResClusterMap.put(track, smallestdrCluster);
            System.out.println("[KFTEM] smallest dr: " + smallestdr);
            System.out.println("[KFTEM] Cluster smallest dr: " + clusterResMap.get(smallestdrCluster));
        }
        return trackMinResClusterMap;
        
    }


}
