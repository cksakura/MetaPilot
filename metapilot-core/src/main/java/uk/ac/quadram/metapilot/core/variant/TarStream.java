package uk.ac.quadram.metapilot.core.variant;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;

/**
 * Minimal sequential tar reader, enough for the UHGP-100 archive.
 * <p>
 * Written here rather than adding a dependency for one archive. Handles what
 * that archive needs: members larger than 8 GiB (the FASTA is 41.6 GB), whose
 * size GNU tar stores in base-256 and POSIX tar in a PAX {@code size} record;
 * long names (GNU {@code L} entries and PAX {@code path}); the ustar name
 * prefix. Members are read in order and cannot be revisited - a gzip stream
 * cannot seek, so reaching the FASTA means reading through the 25 GB TSV
 * before it.
 */
final class TarStream {

	static final class Entry {
		final String name;
		final long size;
		final char type;

		Entry(String name, long size, char type) {
			this.name = name;
			this.size = size;
			this.type = type;
		}

		boolean isFile() {
			return type == '0' || type == '\0' || type == '7';
		}
	}

	private static final int BLOCK = 512;

	private final InputStream in;
	private long remaining;
	private long padding;

	TarStream(InputStream in) {
		this.in = in;
	}

	/** The next member, or null at the end of the archive. Skips what is left of the current one. */
	Entry next() throws IOException {
		skipFully(remaining + padding);
		remaining = 0;
		padding = 0;
		String longName = null;
		Long paxSize = null;
		byte[] h = new byte[BLOCK];
		while (true) {
			if (!readBlock(h)) {
				return null;
			}
			if (allZero(h)) {
				return null;
			}
			char type = (char) h[156];
			long size = parseSize(h);
			long pad = (BLOCK - size % BLOCK) % BLOCK;
			if (type == 'L' || type == 'x' || type == 'g') {
				byte[] data = readFully(size);
				skipFully(pad);
				if (type == 'L') {
					longName = cString(data, 0, data.length);
				} else if (type == 'x') {
					String[] pax = parsePax(data);
					if (pax[0] != null) {
						longName = pax[0];
					}
					if (pax[1] != null) {
						paxSize = Long.parseLong(pax[1]);
					}
				}
				continue;
			}
			String name = longName != null ? longName : headerName(h);
			if (paxSize != null) {
				size = paxSize;
				pad = (BLOCK - size % BLOCK) % BLOCK;
			}
			remaining = size;
			padding = pad;
			return new Entry(name, size, type);
		}
	}

	/** The current member's content; reading past its end returns -1. */
	InputStream content() {
		return new InputStream() {
			@Override
			public int read() throws IOException {
				if (remaining <= 0) {
					return -1;
				}
				int b = in.read();
				if (b < 0) {
					throw new EOFException("Truncated tar member");
				}
				remaining--;
				return b;
			}

			@Override
			public int read(byte[] buf, int off, int len) throws IOException {
				if (remaining <= 0) {
					return -1;
				}
				int n = in.read(buf, off, (int) Math.min(len, remaining));
				if (n < 0) {
					throw new EOFException("Truncated tar member");
				}
				remaining -= n;
				return n;
			}

			@Override
			public void close() {
				// The archive stream is owned by TarStream's caller.
			}
		};
	}

	private static String headerName(byte[] h) {
		String name = cString(h, 0, 100);
		// POSIX "ustar\0" only: old GNU headers ("ustar  \0") keep atime/ctime where the prefix would be.
		boolean posix = h[257] == 'u' && h[258] == 's' && h[259] == 't' && h[260] == 'a' && h[261] == 'r'
				&& h[262] == 0;
		if (posix) {
			String prefix = cString(h, 345, 155);
			if (!prefix.isEmpty()) {
				name = prefix + "/" + name;
			}
		}
		return name;
	}

	/** Octal ASCII, or GNU base-256 when the first byte has its high bit set. */
	static long parseSize(byte[] h) {
		if ((h[124] & 0x80) != 0) {
			long v = h[124] & 0x7f;
			for (int i = 125; i < 136; i++) {
				v = (v << 8) | (h[i] & 0xff);
			}
			return v;
		}
		long v = 0;
		for (int i = 124; i < 136; i++) {
			byte b = h[i];
			if (b == 0 || b == ' ') {
				if (v > 0) {
					break;
				}
				continue;
			}
			if (b < '0' || b > '7') {
				break;
			}
			v = v * 8 + (b - '0');
		}
		return v;
	}

	/** {@code [path, size]} from a PAX extended header; either may be null. */
	private static String[] parsePax(byte[] data) {
		String[] out = new String[2];
		String s = new String(data, StandardCharsets.UTF_8);
		int pos = 0;
		while (pos < s.length()) {
			int space = s.indexOf(' ', pos);
			if (space < 0) {
				break;
			}
			int len;
			try {
				len = Integer.parseInt(s.substring(pos, space));
			} catch (NumberFormatException e) {
				break;
			}
			if (len <= 0 || pos + len > s.length()) {
				break;
			}
			String record = s.substring(space + 1, pos + len - 1);
			int eq = record.indexOf('=');
			if (eq > 0) {
				String key = record.substring(0, eq);
				if (key.equals("path")) {
					out[0] = record.substring(eq + 1);
				} else if (key.equals("size")) {
					out[1] = record.substring(eq + 1);
				}
			}
			pos += len;
		}
		return out;
	}

	private static String cString(byte[] b, int off, int len) {
		int end = off;
		while (end < off + len && end < b.length && b[end] != 0) {
			end++;
		}
		return new String(b, off, end - off, StandardCharsets.UTF_8).trim();
	}

	private static boolean allZero(byte[] b) {
		for (byte x : b) {
			if (x != 0) {
				return false;
			}
		}
		return true;
	}

	private boolean readBlock(byte[] b) throws IOException {
		int off = 0;
		while (off < b.length) {
			int n = in.read(b, off, b.length - off);
			if (n < 0) {
				if (off == 0) {
					return false;
				}
				throw new EOFException("Truncated tar header");
			}
			off += n;
		}
		return true;
	}

	private byte[] readFully(long size) throws IOException {
		if (size > 1 << 24) {
			throw new IOException("Tar metadata entry too large: " + size);
		}
		byte[] b = new byte[(int) size];
		int off = 0;
		while (off < b.length) {
			int n = in.read(b, off, b.length - off);
			if (n < 0) {
				throw new EOFException("Truncated tar metadata entry");
			}
			off += n;
		}
		return b;
	}

	private void skipFully(long n) throws IOException {
		byte[] buf = null;
		while (n > 0) {
			long s = in.skip(n);
			if (s <= 0) {
				if (buf == null) {
					buf = new byte[64 << 10];
				}
				int r = in.read(buf, 0, (int) Math.min(buf.length, n));
				if (r < 0) {
					throw new EOFException("Truncated tar archive");
				}
				s = r;
			}
			n -= s;
		}
	}
}
