package org.hps.record.evio.crawler;

import java.io.File;
import java.io.IOException;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import org.jlab.coda.jevio.EvioException;
import org.jlab.coda.jevio.EvioReader;
import org.lcsim.util.log.LogUtil;

final class EvioFileList extends ArrayList<File> {

    private static final Logger LOGGER = LogUtil.create(EvioFileList.class);

    Map<File, Integer> eventCounts = new HashMap<File, Integer>();

    void computeEventCount(EvioReader reader, final File file) throws IOException, EvioException {
        if (!reader.getPath().equals(file.getPath())) {
            throw new IllegalArgumentException("The EvioReader and file paths do not match.");
        }        
        LOGGER.info("computing event count for " + file.getPath() + " ...");
        int eventCount = reader.getEventCount();
        this.eventCounts.put(file, eventCount);
        LOGGER.info("done computing event count for " + file.getPath());
    }

    int computeTotalEvents() {
        LOGGER.info("computing total events ...");
        int totalEvents = 0;
        for (final File file : this) {
            getEventCount(file);
            totalEvents += this.eventCounts.get(file);
        }
        LOGGER.info("done computing total events");
        return totalEvents;
    }

    File first() {
        return this.get(0);
    }

    int getEventCount(final File file) {
        if (this.eventCounts.get(file) == null) {
            throw new RuntimeException("The event count for " + file.getPath() + " was never computed.");
        }
        return this.eventCounts.get(file);
    }

    void insert(final Connection connection, final int run) throws SQLException {
        LOGGER.info("updating file list ...");
        PreparedStatement filesStatement = null;
        filesStatement = connection
                .prepareStatement("INSERT INTO run_log_files (run, directory, name) VALUES(?, ?, ?)");
        LOGGER.info("inserting files from run " + run + " into database");
        for (final File file : this) {
            LOGGER.info("creating update statement for " + file.getPath());
            filesStatement.setInt(1, run);
            filesStatement.setString(2, file.getParentFile().getPath());
            filesStatement.setString(3, file.getName());
            LOGGER.info("executing statement: " + filesStatement);
            filesStatement.executeUpdate();
        }
        LOGGER.info("run_log_files was updated!");
    }

    File last() {
        return this.get(this.size() - 1);
    }

    void sort() {
        final List<File> fileList = new ArrayList<File>(this);
        Collections.sort(fileList, new EvioFileSequenceComparator());
        this.clear();
        this.addAll(fileList);
    }
}