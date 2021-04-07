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

//Plots for checking new Track MCP matching tools
        plots1D.put(String.format("n_mcps_on_track"), histogramFactory.createHistogram1D(String.format("n_mcps_on_track"), 30, 0, 30));

        plots1D.put(String.format("n_hits"), histogramFactory.createHistogram1D(String.format("n_hits"), 13, 0, 13));

        plots1D.put(String.format("n_goodhits"), histogramFactory.createHistogram1D(String.format("n_goodhits"), 13, 0, 13));

        plots1D.put(String.format("n_badhits"), histogramFactory.createHistogram1D(String.format("n_badhits"), 13, 0, 13));

        plots1D.put(String.format("purity"), histogramFactory.createHistogram1D(String.format("purity"), 100, 0, 1));

        plots1D.put(String.format("layers_hit"), histogramFactory.createHistogram1D(String.format("layers_hit"), 13, 0, 13));

        plots1D.put(String.format("n_mcps_on_layer"), histogramFactory.createHistogram1D(String.format("n_mcps_on_layer"), 20, 0, 20));

        plots1D.put(String.format("n_striphits_on_layer"), histogramFactory.createHistogram1D(String.format("n_striphits_on_layer"), 20, 0, 20));

        plots1D.put(String.format("n_mcps_on_striphit"), histogramFactory.createHistogram1D(String.format("n_mcps_on_striphit"), 10, 0, 10));

        plots2D.put(String.format("n_mcps_per_layer"), histogramFactory.createHistogram2D(String.format("n_mcps_per_layer"), 20, 0, 20, 13, 0, 13));

        plots2D.put(String.format("n_striphits_per_layer"), histogramFactory.createHistogram2D(String.format("n_striphits_per_layer"), 20, 0, 20, 13, 0, 13));

        plots2D.put(String.format("n_mcps_on_layer_striphits"), histogramFactory.createHistogram2D(String.format("n_mcps_on_layer_striphits"), 20, 0, 20, 13, 0, 13));

        plots2D.put(String.format("new_trackMCP_match_trackP_v_mcpP"), histogramFactory.createHistogram2D(String.format("new_trackMCP_match_trackP_v_mcpP"), 1000, -5, 5, 1000, -5, 5));

        plots2D.put(String.format("existing_trackMCP_match_trackP_v_mcpP"), histogramFactory.createHistogram2D(String.format("existing_trackMCP_match_trackP_v_mcpP"), 1000, -5, 5, 1000, -5, 5));

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


        //MCParticleRelations
        
        List<LCRelation> trackToMCParticleRelations    =  new ArrayList<LCRelation>();
        
        //Truth Tracks and Relations
        List<LCRelation> trackToTruthTrackRelations    =  new ArrayList<LCRelation>();
        List<Track>      truthTrackCollection          =  new ArrayList<Track>();
        
        for (Track track : trackCollection) {

            int charge = -1* (int)Math.signum(track.getTrackStates().get(0).getOmega());
            double trackPmag = new BasicHep3Vector(track.getTrackStates().get(0).getMomentum()).magnitude();
            int mcpcharge = 1;

            //Check new truth matching tool
            TrackTruthMatching_new tt = new TrackTruthMatching_new(track, rawtomc);
            MCParticle mcp = tt.getMCParticle();

            if(mcp != null){
                if(mcp.getPDGID() == -11)
                    mcpcharge = 1;
                else
                    mcpcharge = -1;
                plots2D.get("new_trackMCP_match_trackP_v_mcpP").fill(charge*trackPmag,mcpcharge*mcp.getMomentum().magnitude());
                plots1D.get("n_mcps_on_track").fill(tt.getLayerHitsForAllMCPs().size());
                plots1D.get("n_hits").fill(tt.getNHits());
                plots1D.get("purity").fill(tt.getPurity());
                plots1D.get("n_goodhits").fill(tt.getNGoodHits());
                plots1D.get("n_badhits").fill(tt.getNBadHits());
                for(Integer layer : tt.getLayersOnTrack()){
                    plots1D.get("layers_hit").fill(layer);
                    plots1D.get("n_mcps_on_layer").fill(tt.getMCPsOnLayer(layer).size());
                    plots2D.get("n_mcps_per_layer").fill(layer,tt.getMCPsOnLayer(layer).size());
                    plots1D.get("n_striphits_on_layer").fill(tt.getStripHitsOnLayer(layer).size());
                    plots2D.get("n_striphits_per_layer").fill(layer,tt.getStripHitsOnLayer(layer).size());
                    for(RawTrackerHit rawhit : tt.getStripHitsOnLayer(layer)){
                        plots1D.get("n_mcps_on_striphit").fill(tt.getMCPsOnRawTrackerHit(rawhit).size());
                        plots2D.get("n_mcps_on_layer_striphits").fill(layer,tt.getMCPsOnRawTrackerHit(rawhit).size());
                    }
                }
            }

            TrackTruthMatching ttm = new TrackTruthMatching(track, rawtomc, allsimhits, false);
            MCParticle mcp_matt = null;
            if(ttm != null)
                mcp_matt = ttm.getMCParticle();
            if(mcp_matt != null){
                if(mcp_matt.getPDGID() == -11)
                    mcpcharge = 1;
                else
                    mcpcharge = -1;
                plots2D.get("existing_trackMCP_match_trackP_v_mcpP").fill(charge*trackPmag,mcpcharge*mcp_matt.getMomentum().magnitude());
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
}
