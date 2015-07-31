package org.hps.test.it;

import hep.aida.IHistogram1D;

import java.io.File;
import java.io.IOException;
import java.util.List;

import junit.framework.TestCase;

import org.hps.job.JobManager;
import org.hps.recon.ecal.cluster.ClusterUtilities;
import org.hps.test.util.TestFileUrl;
import org.hps.test.util.TestOutputFile;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;
import org.lcsim.util.loop.LCSimLoop;

/**
 * Run readout simulation and full 2015 Engineering Run reconstruction on ECal MC input data and then check histograms of cluster data.
 * 
 * @author Jeremy McCormick, SLAC
 */
public class EcalSimReconTest extends TestCase {
        
    /**
     * Test input file name on lcsim.org site.
     */
    private static final String FILE_NAME = "egsv3-triv2-g4v1_s2d6_HPS-EngRun2015-Nominal-v1_101.slcio";
    
    /**
     * Steering resource file for readout simulation.
     */    
    private static final String READOUT_STEERING = "/org/hps/steering/readout/EngineeringRun2015TrigSingles1.lcsim";
    
    /**
     * Steering resource file for running reconstruction.
     */
    private static final String RECON_STEERING = "/org/hps/steering/recon/EngineeringRun2015FullReconMC.lcsim";
    
    /**
     * Run number for conditions system.
     */
    private static final Integer RUN = 5000;
    
    /**
     * Expected mean of high cluster energy in GeV.
     */
    private static final double HIGH_CLUS_MEAN_E = 1.00;
    
    /**
     * Acceptable delta of high cluster mean energy in GeV.
     */
    private static final double HIGH_CLUS_DELTA_E = 0.05;
    
    /**
     * Expected mean time of primary cluster in nanoseconds.
     */
    private static final double CLUS_MEAN_T = 58.6;   
    
    /**
     * Acceptable delta of cluster time in nanoseconds.
     */    
    private static final double CLUS_DELTA_T = 5;
    
    /**
     * Expected number of triggers.
     */
    private static final int TRIGGERS = 880;
    
    /**
     * Acceptable delta of trigger count.
     */
    private static final int TRIGGER_DELTA = 5;

    
    /**
     * Run the test.
     * 
     * @throws Exception if there is an uncaught exception thrown
     */
    public void testEcalSimRecon() throws Exception {
        
        // Get the input events file.
        File inputFile = TestFileUrl.getInputFile(EcalSimReconTest.class, FILE_NAME);
                
        // Run readout simulation on MC input.
        File readoutFile = new TestOutputFile(EcalSimReconTest.class, "readout.slcio");
        JobManager job = new JobManager();
        job.addVariableDefinition("detector", "HPS-EngRun2015-Nominal-v1");
        job.addVariableDefinition("outputFile", readoutFile.getPath().replace(".slcio", ""));
        job.addVariableDefinition("run", RUN.toString());
        job.addInputFile(inputFile);
        job.setup(READOUT_STEERING);
        job.configure();
        job.run();
        System.out.println("readout ran on " + job.getLCSimLoop().getTotalCountableConsumed() + " events");
        
        // Run the recon on the readout output.
        File reconFile = new TestOutputFile(EcalSimReconTest.class, "recon.slcio");
        job = new JobManager();
        job.addInputFile(readoutFile);
        job.addVariableDefinition("detector", "HPS-EngRun2015-Nominal-v1");
        job.addVariableDefinition("outputFile", reconFile.getPath().replace(".slcio", ""));
        job.addVariableDefinition("run", RUN.toString());
        job.setup(RECON_STEERING);
        job.run();
        System.out.println("recon ran on " + job.getLCSimLoop().getTotalCountableConsumed() + " events");
        
        // Check the recon output.
        LCSimLoop loop = new LCSimLoop();
        loop.setLCIORecordSource(reconFile);
        EcalSimReconCheckDriver checkDriver = new EcalSimReconCheckDriver();
        loop.add(checkDriver);
        loop.loop(-1);
        loop.dispose();
        System.out.println("EcalSimReconCheckDriver ran on " + loop.getTotalCountableConsumed() + " events");
        
        // Check the number of triggers.
        TestCase.assertEquals("Number of triggers out of range.", TRIGGERS, loop.getTotalCountableConsumed(), TRIGGER_DELTA);
    }
    
    /**
     * Driver for checking test output.
     * 
     * @author Jeremy McCormick, SLAC
     */
    static class EcalSimReconCheckDriver extends Driver {
        
        /**
         * Setup AIDA histogramming.
         */
        private final static AIDA aida = AIDA.defaultInstance();
        
        /**
         * Count of clusters per event.
         */
        private final IHistogram1D clusCountH1D = aida.histogram1D("Cluster Count", 20, -0.5, 19.5);
        
        /**
         * Individual cluster energies.
         */
        private final IHistogram1D clusEnergyH1D = aida.histogram1D("Cluster Energy", 100, -0.5, 9.5);
        
        /**
         * Highest cluster energy in event.
         */
        private final IHistogram1D clusHighEnergyH1D = aida.histogram1D("Cluster Highest Energy", 100, -0.5, 9.5);
        
        /**
         * First hit time in highest energy cluster.
         */
        private final IHistogram1D clusTimeH1D = aida.histogram1D("Cluster Time", 500, -0.5, 499.5);        
              
        /**
         * Process events and fill histograms from cluster collection.
         */
        public void process(EventHeader event) {
            List<Cluster> clusters = event.get(Cluster.class, "EcalClustersCorr");
            clusCountH1D.fill(clusters.size());
            for (Cluster cluster : clusters) {
                clusEnergyH1D.fill(cluster.getEnergy());
            }
            
            Cluster highClus = ClusterUtilities.findHighestEnergyCluster(clusters);
            clusHighEnergyH1D.fill(highClus.getEnergy());
            clusTimeH1D.fill(highClus.getCalorimeterHits().get(0).getTime());
        }
        
        /**
         * Save histograms and perform checks on statistics.
         */
        public void endOfData() {
            
            // Print some statistics.
            System.out.println("clusCountH1D mean: " + clusCountH1D.mean());
            System.out.println("clusEnergyH1D mean: " + clusEnergyH1D.mean());
            System.out.println("clusHighEnergyH1D mean: " + clusHighEnergyH1D.mean());
            System.out.println("clusTimeH1D mean: " + clusTimeH1D.mean());
            
            // Save plots to AIDA file.
            try {
                aida.saveAs(new TestOutputFile(EcalSimReconTest.class, "plots.aida"));
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
            
            // Check high cluster energy mean.
            TestCase.assertEquals("High cluster energy was not within range.", HIGH_CLUS_MEAN_E, clusHighEnergyH1D.mean(), HIGH_CLUS_DELTA_E);
            
            // Check high cluster time mean.
            TestCase.assertEquals("High cluster mean time was not within range.", CLUS_MEAN_T, clusTimeH1D.mean(), CLUS_DELTA_T);            
        }
    }    
}
