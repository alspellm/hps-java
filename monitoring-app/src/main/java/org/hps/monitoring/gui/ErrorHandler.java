package org.hps.monitoring.gui;

import java.awt.Component;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;

/**
 * <p>
 * An error handling class which is able to do any of the following,
 * depending on what the caller wants to do with the error.
 * </p>
 * <ul>
 * <li>Print a message</li>
 * <li>Print the stack trace</li>
 * <li>Log to a Logger</li>
 * <li>Show an error dialog</li>
 * <li>Raise an exception</li>
 * <li>Exit the application</li>
 * </ul>
 * </p>
 * It mostly uses the "builder" pattern so that the various handling methods
 * can be easily chained, where appropriate.  Some methods are not available
 * for chaining when it doesn't make sense, e.g. popping up a dialog should
 * happen last so this is not chainable.
 * </p>
 */
class ErrorHandler {

    Logger logger;
    Component component;
    Throwable error;
    String message;
        
    /**
     * Constructor.
     * @param component The GUI component to which this object is assigned.
     * @param logger The logger to which messages will be written.
     */
    ErrorHandler(Component component, Logger logger) {
        this.logger = logger;
        this.component = component;
    }
    
    /**
     * Set the error that occurred.  This should always be called 
     * first in a method chain.
     * @param error The error which is a <code>Throwable</code>.
     * @return This object.
     */
    ErrorHandler setError(Throwable error) {
        this.error = error;
        this.message = error.getMessage();
        return this;
    }
    
    /**
     * Set the error message if it differs from the exception's message.
     * @param message The erro message.
     * @return This object.
     */
    ErrorHandler setMessage(String message) {
        this.message = message;
        return this;
    }
        
    /**
     * Print the full stack trace of the error to System.err.
     * @return This object.
     */
    ErrorHandler printStackTrace() {
        error.printStackTrace();
        return this;
    }
    
    /**
     * Print the error message to System.err.
     * @return This object.
     */
    ErrorHandler printMessage() {        
        System.err.println(message);
        return this;
    }
    
    /**
     * Log the error message to the <code>Logger</code>.
     * @return This object.
     */
    ErrorHandler log() {
        logger.log(Level.SEVERE, message);
        return this;
    }
        
    /**
     * Show an error dialog with the message.
     * @return This object. 
     */
    ErrorHandler showErrorDialog() {
        final Runnable runnable = new Runnable() {
            public void run() {        
                JOptionPane.showMessageDialog(
                        component,
                        error.getMessage(), 
                        "Application Error", 
                        JOptionPane.ERROR_MESSAGE);
            }
        };
        SwingUtilities.invokeLater(runnable);
        return this;
    }
    
    /**
     * Rethrow the error as a <code>RuntimeException</code>.
     * Additional methods cannot be chained to this 
     * as they would not be executed.     
     */
    void raiseException() { 
        throw new RuntimeException(message, error);
    }
    
    /**
     * Exit the application.  This is not chainable for obvious reasons.
     */
    void exit() {
        System.err.println("Fatal error.  Application will exit.");
        System.exit(1);
    }                       
}
