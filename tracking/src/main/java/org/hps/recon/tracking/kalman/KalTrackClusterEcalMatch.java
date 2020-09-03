package org.hps.recon.tracking.kalman;

import hep.aida.IAnalysisFactory;
import hep.aida.IHistogram1D;
import hep.aida.IHistogramFactory;
import hep.aida.ITree;
import hep.aida.ref.rootwriter.RootFileStore;

import org.hps.recon.ecal.cluster.ClusterUtilities;
import org.hps.recon.tracking.TrackUtils;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.List;

import org.lcsim.event.Track;
import org.lcsim.event.TrackState;
import org.lcsim.event.Cluster;
import org.lcsim.event.RelationalTable;

public class KalTrackClusterEcalMatch {


    //Used for plots
    boolean enablePlots;
    private ITree tree;
    private IHistogramFactory histogramFactory;
    private Map<String, IHistogram1D> plots1D;
    private String trackType;

    RelationalTable hitToRotated = null;
    RelationalTable hitToStrips = null;

    public KalTrackClusterEcalMatch(String trackCollectionName) {
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
        plots1D.put(String.format("%s_ElectronTrack-Cluster_dt",trackType), histogramFactory.createHistogram1D(String.format("%s_ElectronTrack-Cluster_dt",trackType), 100, -200, 200));
        plots1D.put(String.format("%s_PositronTrack-Cluster_dt",trackType), histogramFactory.createHistogram1D(String.format("%s_PositronTrack-Cluster_dt",trackType), 100, -200, 200));
        plots1D.put(String.format("%s_PositronTrackTime",trackType), histogramFactory.createHistogram1D(String.format("%s_PositronTrackTime",trackType), 100, -200, 200));
        plots1D.put(String.format("%s_ElectronTrackTime",trackType), histogramFactory.createHistogram1D(String.format("%s_ElectronTrackTime",trackType), 100, -200, 200));
        plots1D.put(String.format("%s_Cluster_Timing_(woffset)",trackType), histogramFactory.createHistogram1D(String.format("%s_Cluster_Timing_(woffset)",trackType), 100, -200, 200));

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
        String trackType = trackCollectionName;
        //have to multiply by -1 for some unknown reason!
        int charge = -1* Charge;
        double trkTime = trackTime;
        double trkxpos;
        double trkypos;
        double trkzpos;
        double dxoffset;

        List<Cluster> clusters = Clusters;

        if(trackType.contains("GBLTracks")) {
            trkxpos = TrackUtils.getTrackStateAtECal(track).getReferencePoint()[1]; 
            trkypos = TrackUtils.getTrackStateAtECal(track).getReferencePoint()[2];
            trkzpos = TrackUtils.getTrackStateAtECal(track).getReferencePoint()[0];
            dxoffset = -4.1;
        }

        else {
            
            TrackState ts_ecal = track.getTrackStates().get(track.getTrackStates().size()-1);
            double[] ts_ecalPos = ts_ecal.getReferencePoint();
            trkxpos = ts_ecalPos[0];
            trkypos = ts_ecalPos[1];
            trkzpos = ts_ecalPos[2];
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
        //double offset = 56;
        double offset = trackClusterTimeOffset;
        double tcut = 4.0;
        double xcut = 10.0;
        double ycut = 10.0;

        for(Cluster cluster : clusters) {
            double clusTime = ClusterUtilities.getSeedHitTime(cluster);
            double dt = clusTime - trackClusterTimeOffset - trkTime;

            double clusxpos = cluster.getPosition()[0];
            double clusypos = cluster.getPosition()[1];
            double cluszpos = cluster.getPosition()[2];
            double dx = clusxpos - trkxpos + dxoffset;
            double dy = clusypos - trkypos;
            double dz = cluszpos - trkzpos;
            System.out.format("%s_dx_%f",trackType,dx);
            System.out.format("%s_dy_%f",trackType,dy);
            System.out.format("%s_dt_%f",trackType,dt);
            double dr = Math.sqrt(Math.pow(clusxpos-trkxpos,2) + Math.pow(clusypos-trkypos,2));

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
                plots1D.get(String.format("%s_PositronTrackTime",trackType)).fill(trkTime);
            }
            else {
                plots1D.get(String.format("%s_ElectronTrackTime",trackType)).fill(trkTime);
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




}
