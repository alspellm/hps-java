package org.hps.readout.ecal.triggerbank;

import java.util.logging.Logger;

import org.hps.recon.ecal.cluster.ClusterDriver;
import org.lcsim.util.log.BasicFormatter;
import org.lcsim.util.log.LogUtil;

/**
 * Abstract class <code>SSPTrigger</code> represents the data output
 * by the SSP for a trigger. Individual implementing classes are expected
 * to handle parsing the trigger bit data as appropriate for their type.
 * 
 * @author Kyle McCarty <mccarty@jlab.org>
 */
public class SSPTrigger {
	// Trigger data.
	protected final int type;
	protected final int time;
	protected final int data;
	
	// Logger to output status messages.
	protected static Logger logger = LogUtil.create(ClusterDriver.class, new BasicFormatter(SSPCluster.class.getSimpleName()));
	
	/**
	 * Instantiates a new <code>SSPTrigger</code> with the indicated
	 * trigger data.
	 * @param type - The type of trigger.
	 * @param time - The time at which the trigger occurred in 4 ns
	 * clock-cycles.
	 * @param data - The trigger bit data.
	 */
	public SSPTrigger(int type, int time, int data) {
		// Log any issues with processing the trigger.
		if(!SSPTriggerFactory.isKnownTriggerType(type)) {
			logger.warning(String.format("Trigger type %d is not recognized.", type));
		}
		
		// Store the trigger data.
		this.type = type;
		this.time = 4 * time;
		this.data = data;
		
		// Note that a trigger was made.
		logger.fine(String.format("A trigger of type %d occurred at time %3d with data %d.",
				type, time, data));
	}
	
	/**
	 * Indicates whether the trigger was reported by the bottom SSP
	 * crate or not.
	 * @return Returns <code>true</code> if the trigger was reported
	 * by the bottom crate and <code>false</code> if it was reported
	 * by the top crate.
	 */
	public boolean isBottom() {
		return false;
	}
	
	/**
	 * Indicates whether the trigger was reported by the top SSP
	 * crate or not.
	 * @return Returns <code>true</code> if the trigger was reported
	 * by the top crate and <code>false</code> if it was reported by
	 * the bottom crate.
	 */
	public boolean isTop() {
		return false;
	}
	
	/**
	 * Gets the raw, unparsed trigger data bank for this trigger.
	 * @return Returns the trigger data bank as an <code>int</code>.
	 */
	public int getData() { return data; }
	
	/**
	 * Gets the type code for the trigger.
	 * @return Returns the trigger type as an <code>int</code>.
	 */
	public int getType() { return type; }
	
	/**
	 * Gets the time at which the trigger occurred.
	 * @return Returns the trigger time as an <code>int</code>.
	 */
	public int getTime() { return time; }
}
