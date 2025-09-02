package uk.ac.quadram.metapilot.core.tools;

import java.io.IOException;
import java.io.InputStream;

public class ResourceLoader {
	
	private static final String RESOURCE_PREFIX = "/resources/";

	public static InputStream load(String fileName) throws IOException {
		InputStream inputStream = ResourceLoader.class.getClassLoader().getResourceAsStream(fileName);

		if (inputStream == null) {
			inputStream = ResourceLoader.class.getResourceAsStream(RESOURCE_PREFIX + fileName);
		}

		if (inputStream == null) {
			// Try alternative resource paths
			inputStream = ResourceLoader.class.getResourceAsStream("/" + fileName);
		}

		if (inputStream == null) {
			// Log the attempted paths for debugging
			System.err.println("ResourceLoader: Attempted to load resource '" + fileName + "' from:");
			System.err.println("  - ClassLoader: " + ResourceLoader.class.getClassLoader());
			System.err.println("  - Class: " + ResourceLoader.class);
			System.err.println("  - Resource prefix: " + RESOURCE_PREFIX);
			throw new IOException("Unable to load resource: " + fileName);
		}

		return inputStream;
	}

	@SuppressWarnings("unused")
	private static boolean isRunningFromJar() {
		return ResourceLoader.class.getResource(RESOURCE_PREFIX) != null;
	}
}
