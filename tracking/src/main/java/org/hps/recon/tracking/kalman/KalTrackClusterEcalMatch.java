package org.hps.recon.tracking.kalman;

import hep.aida.IAnalysisFactory;
import hep.aida.IHistogram1D;
import hep.aida.IHistogramFactory;
import hep.aida.ITree;
import hep.aida.ref.rootwriter.RootFileStore;

import hep.physics.vec.Hep3Vector;
import hep.physics.vec.BasicHep3Vector;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.List;

import org.lcsim.event.Track;
import org.lcsim.event.TrackState;
import org.lcsim.event.Cluster;
import org.hps.recon.tracking.CoordinateTransformations;


public class KalTrackClusterEcalMatch {

    boolean enablePlots;
    private ITree tree;
    private IHistogramFactory histogramFactory;
    private Map<String, IHistogram1D> plots1D;
    private String rootFile = "KalTrackClusterMatch.root";


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

        plots1D.put("Track @ Ecal - Ecal Cluster X Residual",histogramFactory.createHistogram1D("Track @ Ecal - Ecal Cluster X Residual",100, -50, 50));
        plots1D.put("Track @ Ecal - Ecal Cluster Y Residual",histogramFactory.createHistogram1D( "Track @ Ecal - Ecal Cluster Y Residual",100, -50, 50));
        plots1D.put("Track @ Ecal - Ecal Cluster Z Residual",histogramFactory.createHistogram1D( "Track @ Ecal - Ecal Cluster Z Residual",100, -50,50));

        plots1D.put("RK Track @ Ecal - Ecal Cluster X Residual",histogramFactory.createHistogram1D("RK Track @ Ecal - Ecal Cluster X Residual",100, -50, 50));
        plots1D.put("RK Track @ Ecal - Ecal Cluster Y Residual",histogramFactory.createHistogram1D( "RK Track @ Ecal - Ecal Cluster Y Residual",100, -50, 50));
        plots1D.put("RK Track @ Ecal - Ecal Cluster Z Residual",histogramFactory.createHistogram1D( "RK Track @ Ecal - Ecal Cluster Z Residual",100, -50,50));
    }
    public void trackClusterDistance(Track kalmanTrack, List<Cluster> Clusters)  {

        Track track = kalmanTrack;
        List<Cluster> clusters = Clusters;
        TrackState ts_ecal = track.getTrackStates().get(track.getTrackStates().size()-1);
        TrackState ts_ecal_RK = track.getTrackStates().get(track.getTrackStates().size()-2);
        double[] ts_ecalPos = ts_ecal.getReferencePoint();
        Hep3Vector ts_ecalPos_RK = new BasicHep3Vector(ts_ecal_RK.getReferencePoint());
        ts_ecalPos_RK = CoordinateTransformations.transformVectorToDetector(ts_ecalPos_RK);
        System.out.println("LOOK HERE!. ReferencePoint via Kalman:" + ts_ecalPos[0] + " " + ts_ecalPos[1] + " " + ts_ecalPos[2]);
        System.out.println("LOOK HERE!. ReferencePoint via RK:" + ts_ecalPos_RK.x() + " " + ts_ecalPos_RK.y() + " " + ts_ecalPos_RK.z());
        for(Cluster cluster : clusters) {
            double[] clusterPos = cluster.getPosition();
            System.out.println("Cluster Position: ");
            for(int i=0; i < clusterPos.length; i++) {
                System.out.println(clusterPos[i]);
            }

            double dr = Math.sqrt(Math.pow(clusterPos[0]-ts_ecalPos[0],2) + Math.pow(clusterPos[1]-ts_ecalPos[1],2));
            //track cluster diff in x/y via Kalman Extrap
            double xdiff = clusterPos[0]-ts_ecalPos[0];
            double ydiff = clusterPos[1]-ts_ecalPos[1];
            double zdiff = clusterPos[2]-ts_ecalPos[2];
            //via RK Extrap
            double xdiffRK = clusterPos[0]-ts_ecalPos_RK.x();
            double ydiffRK = clusterPos[1]-ts_ecalPos_RK.y();
            double zdiffRK = clusterPos[2]-ts_ecalPos_RK.z();

            if(enablePlots) {
                plots1D.get("Track @ Ecal - Ecal Cluster X Residual").fill(xdiff);
                plots1D.get("Track @ Ecal - Ecal Cluster Y Residual").fill(ydiff);
                plots1D.get("Track @ Ecal - Ecal Cluster Z Residual").fill(zdiff);

                plots1D.get("RK Track @ Ecal - Ecal Cluster X Residual").fill(xdiffRK);
                plots1D.get("RK Track @ Ecal - Ecal Cluster Y Residual").fill(ydiffRK);
                plots1D.get("RK Track @ Ecal - Ecal Cluster Z Residual").fill(zdiffRK);
            }
        }
    }




}
