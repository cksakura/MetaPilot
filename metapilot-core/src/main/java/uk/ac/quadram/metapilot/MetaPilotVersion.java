package uk.ac.quadram.metapilot;

import java.io.InputStream;
import java.util.Optional;
import java.util.jar.Manifest;

/**
 * Resolves the MetaPilot application version for metadata written to result files.
 */
public final class MetaPilotVersion {

	private static final String FALLBACK_VERSION = "1.2.0";

	private MetaPilotVersion() {
	}

	public static String getCurrentVersion() {
		String version = readPackageVersion(MetaPilotVersion.class);
		if (isValid(version)) {
			return version;
		}
		version = readManifestVersion();
		if (isValid(version)) {
			return version;
		}
		return FALLBACK_VERSION;
	}

	private static String readPackageVersion(Class<?> anchor) {
		return Optional.ofNullable(anchor.getPackage())
				.map(Package::getImplementationVersion)
				.orElse(null);
	}

	private static String readManifestVersion() {
		try {
			var url = MetaPilotVersion.class.getClassLoader().getResource("META-INF/MANIFEST.MF");
			if (url == null) {
				return null;
			}
			try (InputStream inputStream = url.openStream()) {
				Manifest manifest = new Manifest(inputStream);
				return manifest.getMainAttributes().getValue("Implementation-Version");
			}
		} catch (Exception ignored) {
			return null;
		}
	}

	private static boolean isValid(String version) {
		return version != null && !version.isBlank();
	}
}
