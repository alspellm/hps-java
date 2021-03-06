package org.hps.readout.util;

import org.hps.readout.ReadoutDriver;

/**
 * Class <code>TriggerTime</code> stores the time at which a trigger
 * occurs in terms of simulation time, and also the triggering driver
 * is tracked.
 */
public class TriggerTime implements Comparable<TriggerTime> {
    /**
     * The simulation time (corrected for time offsets) at which a
     * trigger occurs.
     */
    private final double time;
    /**
     * The driver that produced the trigger.
     */
    private final ReadoutDriver trigger;
    
    /**
     * Trigger type
     */
    private String triggerType = "noSet";
    
    /**
     * Top/bot singles trigger
     */
    private String topBot = "noSet";
    
    /**
     * Instantiates a new <code>TriggerTime</code> object.
     * @param time - The simulation time (corrected for time offsets)
     * at which a trigger occurs.
     * @param trigger - The driver that produced the trigger.
     */
    public TriggerTime(double time, ReadoutDriver trigger) {
        this.time = time;
        this.trigger = trigger;
    }
    
    /**
     * Instantiates a new <code>TriggerTime</code> object.
     * @param time - The simulation time (corrected for time offsets)
     * at which a trigger occurs.
     * @param triggerType - trigger type
     * @param trigger - The driver that produced the trigger.
     */
    public TriggerTime(double time, String triggerType, ReadoutDriver trigger) {
        this.time = time;
        this.trigger = trigger;
        this.triggerType = triggerType;
    }
    
    /**
     * Instantiates a new <code>TriggerTime</code> object.
     * @param time - The simulation time (corrected for time offsets)
     * at which a trigger occurs.
     * @param triggerType - trigger type
     * @param topBot - top/Bot singles trigger
     * @param trigger - The driver that produced the trigger.
     */
    public TriggerTime(double time, String triggerType, String topBot, ReadoutDriver trigger) {
        this.time = time;
        this.triggerType = triggerType;
        this.topBot = topBot;
        this.trigger = trigger;
    }
    
    @Override
    public int compareTo(TriggerTime tt) {
        if(time != tt.time) {
            return Double.compare(time, tt.time);
        } else {
            return 0;
        }
    }
    
    /**
     * Gets the driver class which produced the trigger.
     * @return Returns the class of the trigger driver which created
     * the trigger.
     */
    public ReadoutDriver getTriggeringDriver() {
        return trigger;
    }
    
    /**
     * Gets the simulation time at which the trigger occurred. This
     * value is corrected for the time offsets introduced by any
     * relevant drivers.
     * @return Returns the simulation time of the trigger in units of
     * nanoseconds.
     */
    public double getTriggerTime() {
        return time;
    }
    
    /**
     * Gets trigger type
     */
    public String getTriggerType() {
        return triggerType;
    }
    
    /**
     * Gets top/bot status
     */
    public String getTopBotStat() {
        return topBot;
    }
}
