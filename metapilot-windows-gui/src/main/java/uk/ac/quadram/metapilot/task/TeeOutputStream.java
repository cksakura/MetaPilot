package uk.ac.quadram.metapilot.task;

import java.io.IOException;
import java.io.OutputStream;

/**
 * An OutputStream that writes to two OutputStreams simultaneously.
 * 
 * This class is used to "tee" output - sending it to two destinations at once,
 * similar to the Unix 'tee' command. This is particularly useful for capturing
 * System.out and System.err to both the console and a log file.
 * 
 * Example usage:
 * <pre>
 * FileOutputStream fileOut = new FileOutputStream("output.log");
 * TeeOutputStream teeOut = new TeeOutputStream(System.out, fileOut);
 * System.setOut(new PrintStream(teeOut));
 * </pre>
 * 
 * @author Kai Cheng
 * @version 1.0.0
 */
public class TeeOutputStream extends OutputStream {
    
    private final OutputStream out1;
    private final OutputStream out2;
    
    /**
     * Creates a TeeOutputStream that writes to both output streams.
     * 
     * @param out1 The first output stream
     * @param out2 The second output stream
     */
    public TeeOutputStream(OutputStream out1, OutputStream out2) {
        if (out1 == null || out2 == null) {
            throw new IllegalArgumentException("Output streams cannot be null");
        }
        this.out1 = out1;
        this.out2 = out2;
    }
    
    /**
     * Writes the specified byte to both output streams.
     * 
     * @param b the byte to write
     * @throws IOException if an I/O error occurs
     */
    @Override
    public void write(int b) throws IOException {
        out1.write(b);
        out2.write(b);
    }
    
    /**
     * Writes the specified byte array to both output streams.
     * 
     * @param b the byte array to write
     * @throws IOException if an I/O error occurs
     */
    @Override
    public void write(byte[] b) throws IOException {
        out1.write(b);
        out2.write(b);
    }
    
    /**
     * Writes len bytes from the specified byte array starting at offset off 
     * to both output streams.
     * 
     * @param b the byte array to write
     * @param off the start offset in the data
     * @param len the number of bytes to write
     * @throws IOException if an I/O error occurs
     */
    @Override
    public void write(byte[] b, int off, int len) throws IOException {
        out1.write(b, off, len);
        out2.write(b, off, len);
    }
    
    /**
     * Flushes both output streams, forcing any buffered output bytes to be written out.
     * 
     * @throws IOException if an I/O error occurs
     */
    @Override
    public void flush() throws IOException {
        out1.flush();
        out2.flush();
    }
    
    /**
     * Closes both output streams.
     * 
     * @throws IOException if an I/O error occurs
     */
    @Override
    public void close() throws IOException {
        try {
            out1.close();
        } finally {
            out2.close();
        }
    }
}

