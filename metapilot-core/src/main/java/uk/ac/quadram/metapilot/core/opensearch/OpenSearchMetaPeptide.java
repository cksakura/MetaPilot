package uk.ac.quadram.metapilot.core.opensearch;

import org.dom4j.Element;

import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;

/**
 * Peptide-level aggregation of validated open-search PSMs, in the
 * {@link MetaPeptide} shape the shared MAG post-analysis consumes
 * ({@code DdaMagPostAnalysisService}, protein inference, genome scoring,
 * taxonomy and the report writers).
 * <p>
 * The modified sequence carries the delta mass, so PSMs of the same backbone
 * bearing different modifications aggregate as distinct peptides — which is the
 * point of the workflow. Unmodified forms keep the plain sequence.
 */
public class OpenSearchMetaPeptide extends MetaPeptide {

	private final double deltaMass;
	private final int peakId;
	private final String localizedSequence;

	public OpenSearchMetaPeptide(String sequence, String modSeq, String[] proteins, double[] intensity,
			int[] idenType, double deltaMass, int peakId, String localizedSequence) {
		super(sequence, modSeq, proteins, intensity, idenType);
		this.deltaMass = deltaMass;
		this.peakId = peakId;
		this.localizedSequence = localizedSequence;
	}

	public double getDeltaMass() {
		return deltaMass;
	}

	public int getPeakId() {
		return peakId;
	}

	public String getLocalizedSequence() {
		return localizedSequence;
	}

	/** True when this peptide carries an open-search delta mass. */
	public boolean isModified() {
		return peakId >= 0 && Math.abs(deltaMass) > 0.05;
	}

	@Override
	public Object[] getTableObjects() {
		return new Object[] { getSequence(), getModSeq(), deltaMass, peakId, getScore(), getPEP(),
				getTotalMS2Count(), isDecoy() };
	}

	@Override
	public Element getXmlPepElement() {
		// XML export is not part of this workflow's output; the TSV writers in
		// DdaMagReportExportService are what produce its deliverables.
		return null;
	}
}
