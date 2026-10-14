package uk.ac.quadram.metapilot.core.tools;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;

import org.apache.parquet.io.OutputFile;
import org.apache.parquet.io.PositionOutputStream;

/**
 * A Parquet {@link OutputFile} on the local file system through plain
 * {@code java.nio}. Writing through a Hadoop {@code Path} instead goes via
 * {@code RawLocalFileSystem}, which on Windows sets the new file's permissions
 * with {@code winutils.exe}; without it the writer fails right after creating
 * an empty file.
 */
public final class LocalParquetOutputFile implements OutputFile {

	private static final int BUFFER = 1 << 16;

	private final File file;

	public LocalParquetOutputFile(File file) {
		this.file = file;
	}

	@Override
	public PositionOutputStream create(long blockSizeHint) throws IOException {
		return open(StandardOpenOption.CREATE_NEW);
	}

	@Override
	public PositionOutputStream createOrOverwrite(long blockSizeHint) throws IOException {
		return open(StandardOpenOption.CREATE);
	}

	private PositionOutputStream open(StandardOpenOption mode) throws IOException {
		OutputStream out = new BufferedOutputStream(Files.newOutputStream(file.toPath(), mode,
				StandardOpenOption.WRITE, StandardOpenOption.TRUNCATE_EXISTING), BUFFER);
		return new PositionOutputStream() {
			private long pos;

			@Override
			public long getPos() {
				return pos;
			}

			@Override
			public void write(int b) throws IOException {
				out.write(b);
				pos++;
			}

			@Override
			public void write(byte[] b, int off, int len) throws IOException {
				out.write(b, off, len);
				pos += len;
			}

			@Override
			public void flush() throws IOException {
				out.flush();
			}

			@Override
			public void close() throws IOException {
				out.close();
			}
		};
	}

	@Override
	public boolean supportsBlockSize() {
		return false;
	}

	@Override
	public long defaultBlockSize() {
		return 0;
	}

	@Override
	public String getPath() {
		return file.getAbsolutePath();
	}
}
