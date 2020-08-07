package org.hps.recon.tracking.kalman;

import hep.aida.IAnalysisFactory;
import hep.aida.IHistogram1D;
import hep.aida.IHistogramFactory;
import hep.aida.ITree;
import hep.aida.ref.rootwriter.RootFileStore;

import hep.physics.vec.Hep3Vector;
import hep.physics.vec.BasicHep3Vector;

import org.hps.recon.tracking.CoordinateTransformations;
import org.hps.recon.ecal.cluster.ClusterUtilities;

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
    private String rootFile = "KalTrackClusterMatch.root";

    RelationalTable hitToRotated = null;
    RelationalTable hitToStrips = null;


    public void enablePlots(boolean enablePlots) {
        this.enablePlots = enablePlots;
        if (enablePlots ==true) {
            this.bookHistograms();
        }
    }

    public void saveHistograms() {
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

        plots1D = new HashMap<String, IHistogram1D>();
        tree = IAnalysisFactory.create().createTreeFactory().create();
        histogramFactory = IAnalysisFactory.create().createHistogramFactory(tree);

        //Timing Plots
        plots1D.put("ElectronTrack-Cluster_dt", histogramFactory.createHistogram1D("ElectronTrack-Cluster_dt", 100, -200, 200));
        plots1D.put("PositronTrack-Cluster_dt", histogramFactory.createHistogram1D("PositronTrack-Cluster_dt", 100, -200, 200));
        plots1D.put("PositronTrackTime", histogramFactory.createHistogram1D("PositronTrackTime", 100, -200, 200));
        plots1D.put("ElectronTrackTime", histogramFactory.createHistogram1D("ElectronTrackTime", 100, -200, 200));
        plots1D.put("Cluster_Timing_(woffset)", histogramFactory.createHistogram1D("Cluster_Timing_(woffset)", 100, -200, 200));

        //Kalman Extrapolated Residuals
        plots1D.put("PositronTrack@ECal_ECalCluster_dx",histogramFactory.createHistogram1D("PositronTrack@ECal_ECalCluster_dx",100, -50, 50));
        plots1D.put("PositronTrack@ECal_ECalCluster_dy",histogramFactory.createHistogram1D( "PositronTrack@ECal_ECalCluster_dy",100, -50, 50));
        plots1D.put("PositronTrack@ECal_ECalCluster_dz",histogramFactory.createHistogram1D( "PositronTrack@ECal_ECalCluster_dz",100, -50,50));
        plots1D.put("PositronTrack@ECal_ECalCluster_dr",histogramFactory.createHistogram1D( "PositronTrack@ECal_ECalCluster_dr",100, -50, 150));
        plots1D.put("ElectronTrack@ECal_ECalCluster_dx",histogramFactory.createHistogram1D("ElectronTrack@ECal_ECalCluster_dx",100, -50, 50));
        plots1D.put("ElectronTrack@ECal_ECalCluster_dy",histogramFactory.createHistogram1D( "ElectronTrack@ECal_ECalCluster_dy",100, -50, 50));
        plots1D.put("ElectronTrack@ECal_ECalCluster_dz",histogramFactory.createHistogram1D( "ElectronTrack@ECal_ECalCluster_dz",100, -50,50));
        plots1D.put("ElectronTrack@ECal_ECalCluster_dr",histogramFactory.createHistogram1D( "ElectronTrack@ECal_ECalCluster_dr",100, -50,150));
        
        //RK Extrapolated Residuals
        plots1D.put("RK_PositronTrack@ECal_ECalCluster_dx",histogramFactory.createHistogram1D("RK_PositronTrack@ECal_ECalCluster_dx",100, -50, 50));
        plots1D.put("RK_PositronTrack@ECal_ECalCluster_dy",histogramFactory.createHistogram1D( "RK_PositronTrack@ECal_ECalCluster_dy",100, -50, 50));
        plots1D.put("RK_PositronTrack@ECal_ECalCluster_dz",histogramFactory.createHistogram1D( "RK_PositronTrack@ECal_ECalCluster_dz",100, -50,50));
        plots1D.put("RK_ElectronTrack@ECal_ECalCluster_dx",histogramFactory.createHistogram1D("RK_ElectronTrack@ECal_ECalCluster_dx",100, -50, 50));
        plots1D.put("RK_ElectronTrack@ECal_ECalCluster_dy",histogramFactory.createHistogram1D( "RK_ElectronTrack@ECal_ECalCluster_dy",100, -50, 50));
        plots1D.put("RK_ElectronTrack@ECal_ECalCluster_dz",histogramFactory.createHistogram1D( "RK_ElectronTrack@ECal_ECalCluster_dz",100, -50,50));
    }
    public Cluster getKalmanMatchedCluster(Track kalmanTrackHPS, int Charge, List<Cluster> Clusters, double trackTime, double ecalClusterTimeOffset)  {

        //KalmanTrackHPS
        Track track = kalmanTrackHPS;
        //have to multiply by -1 for some unknown reason!
        int charge = -1* Charge;
        double trkTime = trackTime;

        List<Cluster> clusters = Clusters;

        //Track state at ecal via Kalman extrap
        TrackState ts_ecal = track.getTrackStates().get(track.getTrackStates().size()-1);
        double[] ts_ecalPos = ts_ecal.getReferencePoint();

        //Track state at ecal via RK extrap
        TrackState ts_ecalRK = track.getTrackStates().get(track.getTrackStates().size()-2);
        Hep3Vector ts_ecalPos_RK = new BasicHep3Vector(ts_ecalRK.getReferencePoint());
        ts_ecalPos_RK = CoordinateTransformations.transformVectorToDetector(ts_ecalPos_RK);

        if(enablePlots){
            if (charge > 0) {
                plots1D.get("PositronTrackTime").fill(trkTime);
            }
            else {
                plots1D.get("ElectronTrackTime").fill(trkTime);
            }
        }

        Cluster matchedCluster = null;
        double smallestdt = Double.MAX_VALUE;

        for(Cluster cluster : clusters) {
            double clusTime = ClusterUtilities.getSeedHitTime(cluster);
            double[] clusPos = cluster.getPosition();

            //Track cluster distance 
            double dr = Math.sqrt(Math.pow(clusPos[0]-ts_ecalPos[0],2) + Math.pow(clusPos[1]-ts_ecalPos[1],2));
            //double offset = 56;
            double offset = ecalClusterTimeOffset;

            //track cluster diff in x/y via Kalman Extrap
            double dx = clusPos[0]-ts_ecalPos[0];
            double dy = clusPos[1]-ts_ecalPos[1];
            double dz = clusPos[2]-ts_ecalPos[2];

            //via RK Extrap
            double dxRK = clusPos[0]-ts_ecalPos_RK.x();
            double dyRK = clusPos[1]-ts_ecalPos_RK.y();
            double dzRK = clusPos[2]-ts_ecalPos_RK.z();

            //Track Cluster Time Residual
            double dt = clusTime - offset - trkTime;
            double tcmax = 4.0;
            double tcmin = -4.0;

            double xcmax = 10;
            double xcmin = -10;
            double ycmax = 10;
            double ycmin = -10;


            //Extremely simplified track cluster matching. Cluster that passes
            //position cuts and has closest time is matched. This needs to be
            //updated to a real algorithm.
            if((dt > tcmin && dt < tcmax) && (dx > xcmin && dx < xcmax) && (dy > ycmin && dy < ycmax) ) {
                if(Math.abs(dt) < smallestdt) {
                    smallestdt = Math.abs(dt);
                    matchedCluster = cluster;
                }
            }

            
            if(enablePlots) {
                plots1D.get("Cluster_Timing_(woffset)").fill(clusTime-offset);
                if((dt > tcmin && dt < tcmax) && (dx > xcmin && dx < xcmax) && (dy > ycmin && dy < ycmax) ) {
                    if(charge > 0) {
                        //Time residual plot
                        plots1D.get("PositronTrack-Cluster_dt").fill(dt);

                        //Kalman Extrapolated Residuals
                        plots1D.get("PositronTrack@ECal_ECalCluster_dx").fill(dx);
                        plots1D.get("PositronTrack@ECal_ECalCluster_dy").fill(dy);
                        plots1D.get("PositronTrack@ECal_ECalCluster_dz").fill(dz);
                        plots1D.get("PositronTrack@ECal_ECalCluster_dr").fill(dr);

                        //RK Extrapolated Residuals
                        plots1D.get("RK_PositronTrack@ECal_ECalCluster_dx").fill(dxRK);
                        plots1D.get("RK_PositronTrack@ECal_ECalCluster_dy").fill(dyRK);
                        plots1D.get("RK_PositronTrack@ECal_ECalCluster_dz").fill(dzRK);
                        
                    }
                    else {

                        //Time residual plot
                        plots1D.get("ElectronTrack-Cluster_dt").fill(dt);

                        //Kalman Extrapolated Residuals
                        plots1D.get("ElectronTrack@ECal_ECalCluster_dx").fill(dx);
                        plots1D.get("ElectronTrack@ECal_ECalCluster_dy").fill(dy);
                        plots1D.get("ElectronTrack@ECal_ECalCluster_dz").fill(dz);
                        plots1D.get("ElectronTrack@ECal_ECalCluster_dr").fill(dr);

                        //RK Extrapolated Residuals
                        plots1D.get("RK_ElectronTrack@ECal_ECalCluster_dx").fill(dxRK);
                        plots1D.get("RK_ElectronTrack@ECal_ECalCluster_dy").fill(dyRK);
                        plots1D.get("RK_ElectronTrack@ECal_ECalCluster_dz").fill(dzRK);
                        
                    }
                }
            }
            
        }
        if(matchedCluster == null){
            System.out.println("No matching cluster found for kalman track at ECal");
        }
        return matchedCluster;
    }




}
