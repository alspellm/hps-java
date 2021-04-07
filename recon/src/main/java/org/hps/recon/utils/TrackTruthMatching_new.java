

public class TrackTruthDriver_new {

    private ITree tree;
    private IHistogramFactory histogramFactory;
    private Map<String, IHistogram1D> plots1D;
    private Map<String, IHistogram2D> plots2D;
    boolean enablePlots = false;

    private MCParticle _mcp = null;
    private int _nhits;
    private int _ngoodhits;
    private int _nbadhits;
    private double _purity;
    private Map<MCParticle, int[]> mcParticleStripHits = new HashMap<MCParticle, int[]>();
    private Map<MCParticle, Set<Integer>> mcParticleHitsByLayer = new HashMap<MCParticle, Set<Integer>>();

    protected void enablePlots(boolean input){
        enablePlots = input;
    }

    public TrackTruthMatching_new(Track track, RelationalTable rawtomc){
        doAnalysis(track, rawtomc);
    }

    public void bookHistograms(){
        plots1D = new HashMap<String, IHistogram1D>();
        plots2D = new HashMap<String, IHistogram2D>();
        tree = IAnalysisFactory.create().createTreeFactory().create();
        histogramFactory = IAnalysisFactory.create().createHistogramFactory(tree);

//Plots for checking new Track MCP matching tools
        plots1D.put(String.format("track_max_mcp_hit_multiplicity"), histogramFactory.createHistogram1D(String.format("track_max_mcp_hit_multiplicity"), 100, 0, 100));
        
        plots1D.put(String.format("track_most-next_mcp_hit_multiplicity"), histogramFactory.createHistogram1D(String.format("track_most-next_mcp_hit_multiplicity"), 100, 0, 100));
        
        plots1D.put(String.format("track_number_of_mcps"), histogramFactory.createHistogram1D(String.format("track_number_of_mcps"), 40, 0, 40));
        
        plots2D.put(String.format("track_max_mcp_vs_next_best_mcp_hits"), histogramFactory.createHistogram2D(String.format("track_max_mcp_vs_next_best_mcp_hits"), 40, 0, 40, 40, 0, 40));
        
        plots2D.put(String.format("track_max_mcp_vs_next_best_mcp_hits_same_pdgid"), histogramFactory.createHistogram2D(String.format("track_max_mcp_vs_next_best_mcp_hits_same_pdgid"), 40, 0, 40, 40, 0, 40));
        
        plots2D.put(String.format("track_max_mcp_vs_next_best_mcp_hits_different_pdgid"), histogramFactory.createHistogram2D(String.format("track_max_mcp_vs_next_best_mcp_hits_different_pdgid"), 40, 0, 40, 40, 0, 40));

        plots1D.put(String.format("number_strip_hits_per_sensor_layer"), histogramFactory.createHistogram1D(String.format("number_strip_hits_per_sensor_layer"), 10, 0, 10));
        
        plots1D.put(String.format("number_of_mcps_on_striphits"), histogramFactory.createHistogram1D(String.format("number_of_mcps_on_striphits"), 20, 0, 20));

        plots1D.put(String.format("best_mcp_nStripHits_over_total_nStripHits_on_track"), histogramFactory.createHistogram1D(String.format("best_mcp_nStripHits_over_total_nStripHits_on_track"), 100, 0, 2));

        plots2D.put(String.format("best_mcp_nSensorsHit_v_2nd_best_mcp_nSensorsHit"), histogramFactory.createHistogram2D(String.format("best_mcp_nSensorsHit_v_2nd_best_mcp_nSensorsHit"), 14, 0, 14, 14, 0, 14));

        plots2D.put(String.format("number_mcps_on_si_cluster"), histogramFactory.createHistogram2D(String.format("number_mcps_on_si_cluster"), 13, 0, 13, 10, 0, 10));

        plots2D.put(String.format("new_trackMCP_match_trackP_v_mcpP"), histogramFactory.createHistogram2D(String.format("new_trackMCP_match_trackP_v_mcpP"), 1000, -5, 5, 1000, -5, 5));
        
        plots2D.put(String.format("existing_trackMCP_match_trackP_v_mcpP"), histogramFactory.createHistogram2D(String.format("existing_trackMCP_match_trackP_v_mcpP"), 1000, -5, 5, 1000, -5, 5));

    }

    public MCParticle getMCPsOnTrack(Track track, RelationalTable rawtomc){

        double trackPmag = new BasicHep3Vector(track.getTrackStates().get(0).getMomentum()).magnitude();
        int _nhits = track.getTrackerHits().size();

        //loop over tracker hits
        for (TrackerHit hit: track.getTrackerHits()) {
            //Get layer by layer truth information
            Map<Integer, Map<RawTrackerHit, List<MCParticle>>> hitTruth_byLayer = getMCParticlesOnTrackerHit(hit, rawtomc);
            //loop over each layer that this trackerhit exists on
            for(Map.Entry<Integer, Map<RawTrackerHit, List<MCParticle>>> entry : hitTruth_byLayer.entrySet()){
                //loop over each strip hit that was clustered into this layer
                //hit
                for(Map.Entry<RawTrackerHit, List<MCParticle>> subentry : entry.getValue().entrySet()){
                    //loop over each unique MCParticle that contributed to this
                    //strip hit
                    for(MCParticle particle : subentry.getValue()){
                        //Only count a MCParticle hit once per layer
                        Set<Integer> hitsOnLayer = new HashSet<Integer>();
                        hitsOnLayer.add(entry.getKey());
                        if(!mcParticleHitsByLayer.containsKey(particle)){
                            mcParticleHitsByLayer.put(particle, hitsOnLayer);
                        }
                        else{
                            Set<Integer> tmp = mcParticleHitsByLayer.get(particle);
                            tmp.add(entry.getKey());
                            mcParticleHitsByLayer.put(particle, tmp);
                        }
                        if(!mcParticleStripHits.containsKey(particle)){
                            mcParticleStripHits.put(particle, new int[1]);
                            mcParticleStripHits.get(particle)[0] = 0;
                        }
                        mcParticleStripHits.get(particle)[0]++;
                    }
                }
            }
        }
    }

    //match track to MCP by matching to MCP that leaves hits on the most layers
    private MCParticle matchTrackToMCP(Map<MCParticle, Set<Integer>> mcParticleHitsByLayer){

        MCParticle bestMCP = null;
        int bestNhits = 0;
        int totalNLayersHit = 0;
        int allMCPLayerHits = 0;
        for(Map.Entry<MCParticle, Set<Integer>> entry : mcParticleHitsByLayer.entrySet()){
            int nhits = entry.getValue().size();
            if(nhits > bestNhits){
                bestNhits = nhits;
                bestMCP = entry.getKey();
            }
        }

        return bestMCP;
        
    }


    public Map<Integer, Map<RawTrackerHit, List<MCParticle>>> getMCParticlesOnTrackerHit(TrackerHit hit, RelationalTable rawtomc){

        //Get all 1d Strip Hits on TrackerHit
        List<RawTrackerHit> rawhits = hit.getRawHits();
        System.out.println("3d hit made of " + rawhits.size() + " RawTrackerHits");
        //GBL TrackerHits are 3d clusters of 2d hits from axial and stereo pair
        //Separate these 3d clusters into their independent layer 2d constituents
        Set<Integer> layers = new HashSet<Integer>();
        //Define which layers this TrackerHit is composed of
        for(RawTrackerHit rawhit : rawhits){
            int layer = rawhit.getLayerNumber();
            layers.add(layer);
        }

        //Map StripHits to their respective layers
        Map<Integer, Map<RawTrackerHit, List<MCParticle>>> rawHitsMCPMap_byLayer = new HashMap<Integer,Map<RawTrackerHit, List<MCParticle>>>();

        //Loop over Striphits, layer by layer
        for(Integer layer : layers){
            Map<RawTrackerHit, List<MCParticle>> rawHitsMCPMap = new HashMap<RawTrackerHit, List<MCParticle>>();
            System.out.println("Checking RawTrackerHits on layer " + layer);
            int nlayerRawHits = 0;
            int nlayerMCPs = 0;
            for(RawTrackerHit rawhit : rawhits){
                if(rawhit.getLayerNumber() != layer)
                    continue;
                nlayerRawHits = nlayerRawHits + 1;
                //Get list of unique MCPs that make up this Strip hit
                List<MCParticle> rawhitMCPs = getRawHitMCPs(rawhit, rawtomc);
                nlayerMCPs = nlayerMCPs + rawhitMCPs.size();
                rawHitsMCPMap.put(rawhit,rawhitMCPs);
            }
            System.out.println("Total N RawTrackerHits on layer " + layer + " is " +  rawHitsMCPMap.size());
            System.out.println("Total N unique MCPs on layer " + layer + " is " + nlayerMCPs);
            plots1D.get("number_strip_hits_per_sensor_layer").fill(nlayerRawHits);

            rawHitsMCPMap_byLayer.put(layer, rawHitsMCPMap);
            plots2D.get("number_mcps_on_si_cluster").fill(layer,nlayerMCPs);
        }

        return rawHitsMCPMap_byLayer;
    }

    public List<MCParticle> getRawHitMCPs(RawTrackerHit rawhit, RelationalTable rawtomc){

        Set<MCParticle> mcps = new HashSet<MCParticle>();
        List<MCParticle> mcpList = new ArrayList<MCParticle>();
        Set<SimTrackerHit> simhits = rawtomc.allFrom(rawhit);
        System.out.println("Number of simhits on rawtrackerhit: " + simhits.size());
        //loop over simhits on rawhit
        for(SimTrackerHit simhit : simhits){
            //get mcp that left simhit
            MCParticle particle = simhit.getMCParticle();
            int pdgid = particle.getPDGID();
            //only add electron or positron mcps
            if(Math.abs(pdgid) != 11)
                continue;
            if(particle.getOriginZ() > 0)
                continue;
            System.out.println("Found MCP PDGID " + pdgid + " with energy " + particle.getEnergy() + " on RawTrackerHit");
            mcps.add(particle);
        }

        //Convert set to list, so that list contains only 1 entry per MCP
        mcpList.addAll(mcps);
        System.out.println("N MCPs on RawTrackerHit: " + mcps.size());
        plots1D.get("number_of_mcps_on_striphits").fill(mcps.size());

        return mcpList;
    }


    public void plotTrackMCPMultiplicity(Map<MCParticle, int[]> mcParticleStripHits){


        // Look for the MC particle that occurs the most of the track
        int maxValue = 0;
        MCParticle maxMCP = null;
        for(Map.Entry<MCParticle, int[]> entry : mcParticleStripHits.entrySet()){
            if(maxValue < entry.getValue()[0]){
                maxMCP = entry.getKey();
                maxValue = entry.getValue()[0];
            }
        }

        int secondbestHits = 0;
        MCParticle secondbestMCP = null;
        for(Map.Entry<MCParticle, int[]> entry : mcParticleStripHits.entrySet()){
            if(entry.getKey() == maxMCP)
                continue;
            if(secondbestHits < entry.getValue()[0]){
                secondbestMCP = entry.getKey();
                secondbestHits = entry.getValue()[0];
            }
        }

        plots1D.get("track_max_mcp_hit_multiplicity").fill(maxValue);
        plots1D.get("track_most-next_mcp_hit_multiplicity").fill(maxValue - secondbestHits);
        plots1D.get("track_number_of_mcps").fill(mcParticleStripHits.size());
        plots2D.get("track_max_mcp_vs_next_best_mcp_hits").fill(maxValue, secondbestHits);

        if(secondbestMCP != null && maxMCP != null){
            if(maxMCP.getPDGID() == secondbestMCP.getPDGID())
                plots2D.get("track_max_mcp_vs_next_best_mcp_hits_same_pdgid").fill(maxValue, secondbestHits);
            else{ 
                plots2D.get("track_max_mcp_vs_next_best_mcp_hits_different_pdgid").fill(maxValue, secondbestHits);
                if(secondbestHits/maxValue > 0.6){
                    double secondP = secondbestMCP.getMomentum().magnitude();
                    double firstP = maxMCP.getMomentum().magnitude();
                    System.out.println("Best MCP pdgid: " + maxMCP.getPDGID());
                    System.out.println("Best MCP momentum: " + firstP);
                    System.out.println("Second Best MCP pdgid: " + secondbestMCP.getPDGID());
                    System.out.println("Second Best MCP momentum: " + secondP);
                }
            }
        }

        MCParticle bestMCP = null;
        MCParticle secondbestMCP = null;
        int maxHits = 0;
        int secondbestHits = 0;
        int totalHits = 0;

        //find MCP that hits the most layers
        for(Map.Entry<MCParticle, int[]> entry : mcParticleStripHits.entrySet()){
            totalHits = totalHits + entry.getValue()[0];
            if(entry.getValue()[0] > maxHits){
                bestMCP = entry.getKey();
                maxHits = entry.getValue()[0];
            }
            System.out.println("MCP pdgid " + entry.getKey().getPDGID() + " and momentum " + entry.getKey().getMomentum().magnitude() + " leaves "  + entry.getValue()[0] + " strip hits on Track: ");
            System.out.println(" and it hits " + mcParticleHitsByLayer.get(entry.getKey()).size() + " different sensors");
        }

        if(enablePlots){

            for(Map.Entry<MCParticle, int[]> entry : mcParticleStripHits.entrySet()){
                if(entry.getKey() == bestMCP)
                    continue;
                if(secondbestHits < entry.getValue()[0]){
                    secondbestMCP = entry.getKey();
                    secondbestHits = entry.getValue()[0];
                }
            }

            plots1D.get("best_mcp_nStripHits_over_total_nStripHits_on_track").fill((double)maxHits/(double)totalHits);
            if(secondbestMCP != null)
                plots2D.get("best_mcp_nSensorsHit_v_2nd_best_mcp_nSensorsHit").fill(mcParticleHitsByLayer.get(bestMCP).size(),mcParticleHitsByLayer.get(secondbestMCP).size());
            else
                plots2D.get("best_mcp_nSensorsHit_v_2nd_best_mcp_nSensorsHit").fill(mcParticleHitsByLayer.get(bestMCP).size(),0);
        }
    }
}



