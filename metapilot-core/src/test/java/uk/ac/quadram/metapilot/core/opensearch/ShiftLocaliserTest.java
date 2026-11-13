package uk.ac.quadram.metapilot.core.opensearch;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import uk.ac.quadram.metapilot.core.variant.MzmlScanReader;

class ShiftLocaliserTest {

	private static final String PEPTIDE = "LGEYGFQNSLIVR";
	private static final int S = PEPTIDE.indexOf('S');
	private static final double CH2 = 14.01565;

	@TempDir
	File tmp;

	@Test
	@DisplayName("A CH2 shift carried by the fragments of the serine is localised there and explained (S->T)")
	void localisesAndExplains() {
		ShiftLocaliser.Result r = new ShiftLocaliser(20).localise(PEPTIDE, 2, CH2, peaks(PEPTIDE, S, CH2, 1));
		assertEquals(ShiftLocaliser.State.LOCALISED_EXPLAINED, r.state);
		assertEquals(S, r.site);
		assertTrue(r.siteProbability > 0.9, "site probability " + r.siteProbability);
		assertTrue(r.support > 5, "support " + r.support);
		assertTrue(r.explanations.contains("subst S->T"), r.explanations.toString());
	}

	@Test
	@DisplayName("A shift no fragment carries stays ambiguous, with nothing claimed")
	void unshiftedSpectrumIsAmbiguous() {
		ShiftLocaliser.Result r = new ShiftLocaliser(20).localise(PEPTIDE, 2, CH2, peaks(PEPTIDE, -1, CH2, 2));
		assertEquals(ShiftLocaliser.State.UNSUPPORTED, r.state);
		assertTrue(r.support <= 0, "support " + r.support);
		assertTrue(r.explanations.isEmpty());
	}

	@Test
	@DisplayName("A confidently localised mass that nothing explains at its site is localised but unexplained")
	void unknownMassAtSite() {
		double odd = 7.7777;
		int site = PEPTIDE.indexOf('Y');
		ShiftLocaliser.Result r = new ShiftLocaliser(20).localise(PEPTIDE, 2, odd, peaks(PEPTIDE, site, odd, 3));
		assertEquals(ShiftLocaliser.State.LOCALISED_UNEXPLAINED, r.state);
		assertEquals(site, r.site);
		assertTrue(r.explanations.isEmpty());
	}

	private static final double DEAMIDATION = 0.98402;
	private static final double PYRO_GLU = -17.02655;

	@Test
	@DisplayName("Accepted-shift check: no fragment support, consistent explained residue, contradicted residues")
	void acceptedShifts() {
		ShiftLocaliser loc = new ShiftLocaliser(20);
		// Deamidation carried by the fragments of the N: kept, on the N.
		int n = PEPTIDE.indexOf('N');
		ShiftLocaliser.Result onN = loc.checkAcceptedShift(PEPTIDE, 2, DEAMIDATION, peaks(PEPTIDE, n, DEAMIDATION, 11));
		assertEquals(ShiftLocaliser.State.LOCALISED_EXPLAINED, onN.state);
		assertEquals(n, onN.site);
		assertTrue(onN.explanations.stream().anyMatch(e -> e.contains("Deamidat")), onN.explanations.toString());
		assertEquals(DEAMIDATION, onN.explainedMass, 1e-4);
		// No shifted fragments: unsupported, so exported unmodified.
		assertEquals(ShiftLocaliser.State.UNSUPPORTED,
				loc.checkAcceptedShift(PEPTIDE, 2, DEAMIDATION, peaks(PEPTIDE, -1, DEAMIDATION, 12)).state);
		// +0.984 carried by the fragments of a peptide without N or Q: no explained residue, so no explanation.
		String noNQ = "LGEYGFAKALIVR";
		ShiftLocaliser.Result contradicted = loc.checkAcceptedShift(noNQ, 2, DEAMIDATION,
				peaks(noNQ, 6, DEAMIDATION, 13));
		assertEquals(ShiftLocaliser.State.LOCALISED_UNEXPLAINED, contradicted.state);
		assertTrue(contradicted.explanations.isEmpty());
	}

	private static final double A_TO_T = 30.010565;
	private static final double OXIDATION = 15.994915;

	@Test
	@DisplayName("UNIMOD's amino acid substitution entries are no modification: a substitution alone is a "
			+ "substitution, and only a real modification on the same residue makes it undecidable")
	void explanationClasses() {
		ShiftLocaliser loc = new ShiftLocaliser(20);
		// +30.0106 on A: A->T, and nothing else there but UNIMOD's own "Ala->Thr".
		String withA = "LGEYGFANSLIVR";
		int a = withA.indexOf('A');
		ShiftLocaliser.Result subst = loc.localise(withA, 2, A_TO_T, peaks(withA, a, A_TO_T, 21));
		assertEquals(ShiftLocaliser.State.LOCALISED_EXPLAINED, subst.state);
		assertEquals(List.of("subst A->T"), subst.explanations);
		assertEquals(ShiftLocaliser.ExplanationClass.SUBSTITUTION, subst.explanationClass());
		assertEquals('T', subst.substitutedResidue());
		assertEquals("substitution", subst.explanationClassText());
		// CH2 on S: S->T, and methylation of S is a UNIMOD modification, so the spectrum cannot decide.
		ShiftLocaliser.Result both = loc.localise(PEPTIDE, 2, CH2, peaks(PEPTIDE, S, CH2, 1));
		assertTrue(both.explanations.contains("subst S->T") && both.explanations.contains("unimod Methyl"),
				both.explanations.toString());
		assertFalse(both.explanations.contains("unimod Ser->Thr"), both.explanations.toString());
		assertEquals(ShiftLocaliser.ExplanationClass.MODIFICATION_OR_SUBSTITUTION, both.explanationClass());
		// Oxidation on M: no residue is 15.995 heavier than M.
		String withM = "LGEYGFQNMLIVR";
		int m = withM.indexOf('M');
		ShiftLocaliser.Result mod = loc.localise(withM, 2, OXIDATION, peaks(withM, m, OXIDATION, 22));
		assertEquals(ShiftLocaliser.ExplanationClass.MODIFICATION, mod.explanationClass(), mod.explanations.toString());
		assertEquals(0, mod.substitutedResidue());
		assertNull(loc.localise(PEPTIDE, 2, CH2, peaks(PEPTIDE, -1, CH2, 2)).explanationClass());
	}

	@Test
	@DisplayName("An accepted-shift winner explained on a substitution site is set aside like an off-peak one; "
			+ "one on a residue no substitution fits stays identified")
	void acceptedShiftOnSubstitutionSiteIsSetAside() throws IOException {
		String withM = "LGEYGFQNMLIVR";
		List<OpenSearchPsm> open = new ArrayList<>();
		open.add(openPsm(1, PEPTIDE, true, CH2)); // CH2 on S: S->T fits
		open.add(openPsm(2, withM, true, OXIDATION)); // oxidation on M: no substitution fits
		open.add(openPsm(3, PEPTIDE, true, CH2)); // unshifted spectrum: unsupported
		for (OpenSearchPsm p : open) {
			p.setPeakId(5);
		}
		File mzml = new File(tmp, "run1.mzML");
		writeMzml(mzml, Map.of(1, peaks(PEPTIDE, S, CH2, 4), 2, peaks(withM, withM.indexOf('M'), OXIDATION, 22), 3,
				peaks(PEPTIDE, -1, CH2, 6)));
		List<JointCandidate> candidates = JointCandidate.merge(open, Collections.emptyList(), new int[2]);
		for (JointCandidate c : candidates) {
			c.group = JointCandidate.Group.OPEN_ACCEPTED_SHIFT;
		}
		JointIdentification.Config cfg = new JointIdentification.Config();
		cfg.shiftLocaliser = new ShiftLocaliser(20);
		cfg.localisationSpectra.put("run1", mzml);
		List<String> messages = new ArrayList<>();
		JointIdentification.localiseWinners(candidates, new boolean[] { true, true, true }, cfg, messages::add,
				messages::add);
		assertEquals(JointCandidate.Group.OPEN_SUBSTITUTION_SITE, candidates.get(0).getGroup());
		assertEquals("OPEN_SUBSTITUTION_SITE", candidates.get(0).openPrimary().getIdentificationGroup());
		assertEquals(JointCandidate.Group.OPEN_ACCEPTED_SHIFT, candidates.get(1).getGroup());
		assertEquals(ShiftLocaliser.State.LOCALISED_EXPLAINED, candidates.get(1).getLocalisation().state);
		assertEquals(JointCandidate.Group.OPEN_ACCEPTED_SHIFT, candidates.get(2).getGroup(),
				"unplaced, and UNIMOD allows a methyl on this peptide: identified");
		assertEquals(ShiftLocaliser.State.ASSIGNED, candidates.get(2).getLocalisation().state);
		assertTrue(messages.stream().anyMatch(s -> s.contains("1 accepted-shift winners on a substitution site")),
				messages.toString());
	}

	@Test
	@DisplayName("An unplaced accepted shift takes the UNIMOD modification of its mass, on a common site before a "
			+ "rare one; a mass UNIMOD has on no residue of the peptide is not accounted for")
	void unplacedShiftsAreAssignedFromUnimod() {
		ShiftLocaliser loc = new ShiftLocaliser(20);
		// Oxidation: UNIMOD lists L, G, E, ... as rare sites and M as a common one. The M takes it, not residue 1.
		String withM = "LGEYGFQNMLIVR";
		ShiftLocaliser.Result onM = loc.assign(withM, OXIDATION, -1.5);
		assertEquals(ShiftLocaliser.State.ASSIGNED, onM.state);
		assertEquals(withM.indexOf('M'), onM.site);
		assertTrue(onM.explanations.contains("unimod Oxidation"), onM.explanations.toString());
		assertEquals(OXIDATION, onM.explainedMass, 1e-6);
		assertTrue(Double.isNaN(onM.siteProbability), "the site is from the table, not the spectrum");
		assertTrue(onM.isExplained());
		assertEquals(ShiftLocaliser.ExplanationClass.MODIFICATION, onM.explanationClass());
		// Without a common site, the first rare one.
		ShiftLocaliser.Result rare = loc.assign(PEPTIDE, OXIDATION, 0.0);
		assertEquals(ShiftLocaliser.State.ASSIGNED, rare.state);
		assertEquals(0, rare.site);
		// The +4 Da band: no UNIMOD entry.
		assertNull(loc.assign(withM, 4.0531, 0.0));
		assertNull(loc.assign(withM, Double.NaN, 0.0));
	}

	@Test
	@DisplayName("Accepted-shift winners by PSM: unplaced with a UNIMOD site stays identified as that modified "
			+ "peptide, unplaced without one and contradicted ones are set aside")
	void acceptedShiftsFollowTheUnimodRule() throws IOException {
		String withM = "LGEYGFQNMLIVR";
		String noNQ = "LGEYGFAKALIVR";
		double band = 4.0531;
		List<OpenSearchPsm> open = new ArrayList<>();
		open.add(openPsm(1, withM, true, OXIDATION)); // unshifted spectrum, oxidation: assigned to the M
		open.add(openPsm(2, withM, true, band)); // unshifted spectrum, no UNIMOD mass: set aside
		open.add(openPsm(3, noNQ, true, DEAMIDATION)); // +0.984 carried by a residue that cannot deamidate
		for (OpenSearchPsm p : open) {
			p.setPeakId(5);
			p.setCalcNeutralPepMass(1500.0);
		}
		File mzml = new File(tmp, "run1.mzML");
		writeMzml(mzml, Map.of(1, peaks(withM, -1, OXIDATION, 31), 2, peaks(withM, -1, band, 32), 3,
				peaks(noNQ, 6, DEAMIDATION, 13)));
		List<JointCandidate> candidates = JointCandidate.merge(open, Collections.emptyList(), new int[2]);
		for (JointCandidate c : candidates) {
			c.group = JointCandidate.Group.OPEN_ACCEPTED_SHIFT;
		}
		JointIdentification.Config cfg = new JointIdentification.Config();
		cfg.shiftLocaliser = new ShiftLocaliser(20);
		cfg.localisationSpectra.put("run1", mzml);
		List<String> messages = new ArrayList<>();
		JointIdentification.localiseWinners(candidates, new boolean[] { true, true, true }, cfg, messages::add,
				messages::add);

		OpenSearchPsm oxidised = candidates.get(0).openPrimary();
		assertEquals(JointCandidate.Group.OPEN_ACCEPTED_SHIFT, candidates.get(0).getGroup());
		assertEquals(ShiftLocaliser.State.ASSIGNED, oxidised.getShiftLocalisation().state);
		assertEquals(withM + "[+15.9949]", oxidised.getModificationFamilyKey());
		assertEquals(1500.0 + OXIDATION, oxidised.getFlashLfqMonoisotopicMass(), 1e-6,
				"FlashLFQ looks for the oxidised peptide, the precursor that was fragmented");
		assertFalse(oxidised.isShiftFoldedUnmodified());
		assertEquals(withM.indexOf('M'), oxidised.getLocalisedShiftSite());

		assertEquals(JointCandidate.Group.OPEN_UNEXPLAINED, candidates.get(1).getGroup());
		assertEquals(ShiftLocaliser.State.UNSUPPORTED, candidates.get(1).getLocalisation().state);
		assertEquals(JointCandidate.Group.OPEN_UNEXPLAINED, candidates.get(2).getGroup());
		assertEquals(ShiftLocaliser.State.LOCALISED_UNEXPLAINED, candidates.get(2).getLocalisation().state);
		assertFalse(JointIdentification.identifies(candidates.get(1).getGroup(), cfg));
	}

	private static OpenSearchPsm called(String group, List<String> explanations, String oldVerdict) {
		OpenSearchPsm p = openPsm(1, "SAMPLEK", true, A_TO_T);
		p.setIdentificationGroup(group);
		p.setShiftLocalisation(new ShiftLocaliser.Result(8.0, 1, 0.97, ShiftLocaliser.State.LOCALISED_EXPLAINED,
				explanations, A_TO_T));
		if (oldVerdict != null) {
			p.setSubstitutionVerdict(oldVerdict);
			p.setSubstitutionCandidate("TAMPLEK");
			p.setDeltaMassReassignment("SUBSTITUTION:S1T");
		}
		return p;
	}

	@Test
	@DisplayName("Substitution verdicts follow the identification: set on set-aside PSMs a substitution alone "
			+ "explains, removed everywhere else, untouched without an identification group")
	void verdictsFollowTheIdentification() {
		OpenSearchPsm substitution = called("OPEN_SUBSTITUTION_SITE", List.of("subst A->T"), null);
		OpenSearchPsm undecidable = called("OPEN_SUBSTITUTION_SITE", List.of("subst A->T", "unimod Something"),
				"SUBSTITUTION_NOT_IN_DB");
		OpenSearchPsm identified = called("OPEN_ACCEPTED_SHIFT", List.of("unimod Something"), "SUBSTITUTION_IN_DB");
		OpenSearchPsm otherMode = called(null, List.of("subst A->T"), "SUBSTITUTION_IN_DB");
		SubstitutionReassigner.Alignment a = SubstitutionReassigner.alignWithIdentification(
				Arrays.asList(substitution, undecidable, identified, otherMode), null);
		assertEquals(1, a.called);
		assertEquals(0, a.inDatabase);
		assertEquals(2, a.cleared);
		assertEquals("SUBSTITUTION_NOT_IN_DB", substitution.getSubstitutionVerdict());
		assertEquals("STMPLEK", substitution.getSubstitutionCandidate());
		assertEquals("SUBSTITUTION:A2T", substitution.getDeltaMassReassignment());
		assertNull(undecidable.getSubstitutionVerdict());
		assertNull(undecidable.getSubstitutionCandidate());
		assertNull(identified.getSubstitutionVerdict());
		assertEquals("SUBSTITUTION_IN_DB", otherMode.getSubstitutionVerdict());
		assertEquals("TAMPLEK", otherMode.getSubstitutionCandidate());
	}

	@Test
	@DisplayName("An N-terminal pyroGlu without b1 and y(n-1) shares its probability with residue 2, and still stands")
	void pyroGluWithoutB1() {
		String q = "QGEYGFANSLIVR";
		MzmlScanReader.Peaks all = peaks(q, 0, PYRO_GLU, 14);
		// The only ions that tell residue 1 from residue 2: b1 and y(n-1).
		double b1 = SubstitutionMassTable.residueMass('Q') + PYRO_GLU + ShiftLocaliser.PROTON;
		double yN1 = ShiftLocaliserTest.mass(q.substring(1)) + ShiftLocaliser.H2O + ShiftLocaliser.PROTON;
		List<Integer> keep = new ArrayList<>();
		for (int i = 0; i < all.mz.length; i++) {
			if (Math.abs(all.mz[i] - b1) > 0.01 && Math.abs(all.mz[i] - yN1) > 0.01) {
				keep.add(i);
			}
		}
		double[] mz = new double[keep.size()];
		float[] in = new float[keep.size()];
		for (int i = 0; i < mz.length; i++) {
			mz[i] = all.mz[keep.get(i)];
			in[i] = all.intensity[keep.get(i)];
		}
		ShiftLocaliser loc = new ShiftLocaliser(20);
		ShiftLocaliser.Result r = loc.checkAcceptedShift(q, 2, PYRO_GLU, new MzmlScanReader.Peaks(mz, in));
		assertEquals(ShiftLocaliser.State.LOCALISED_EXPLAINED, r.state, "residues 1 and 2 cannot be told apart");
		assertEquals(0, r.site);
		assertTrue(r.siteProbability < 0.75, "a site-probability threshold would have rejected it: "
				+ r.siteProbability);
	}

	@Test
	@DisplayName("Ion coverage of a site and the binomial evidence")
	void ionsAndSignificance() {
		// n = 5: b2 covers residues 0-1, y2 residues 3-4
		assertTrue(ShiftLocaliser.carriesSite(true, 2, 1, 5));
		assertFalse(ShiftLocaliser.carriesSite(true, 2, 2, 5));
		assertTrue(ShiftLocaliser.carriesSite(false, 2, 3, 5));
		assertFalse(ShiftLocaliser.carriesSite(false, 2, 2, 5));
		assertEquals(0.0, ShiftLocaliser.significance(0, 10, 0.1), 0.0);
		assertEquals(-Math.log10(1 - Math.pow(0.9, 10)), ShiftLocaliser.significance(1, 10, 0.1), 1e-9);
		assertTrue(ShiftLocaliser.significance(8, 10, 0.1) > ShiftLocaliser.significance(4, 10, 0.1));
		assertTrue(ShiftLocaliser.significance(4, 10, 0.05) > ShiftLocaliser.significance(4, 10, 0.2),
				"a denser spectrum makes the same matches weaker evidence");
	}

	@Test
	@DisplayName("Explained winners on a substitution site, target or decoy, are reported but not identified; others untouched")
	void localiseWinnersMovesExplainedWinners() throws IOException {
		String decoy = "LVIGSNQFGYEGR";
		List<OpenSearchPsm> open = new ArrayList<>();
		open.add(openPsm(1, PEPTIDE, true, CH2)); // shifted spectrum: explained
		open.add(openPsm(2, decoy, false, CH2)); // decoy, shifted spectrum: explained too (label-blind)
		open.add(openPsm(3, PEPTIDE, true, CH2)); // unshifted spectrum: ambiguous
		open.add(openPsm(4, PEPTIDE, true, CH2)); // not a winner: never localised
		File mzml = new File(tmp, "run1.mzML");
		writeMzml(mzml, Map.of(1, peaks(PEPTIDE, S, CH2, 4), 2, peaks(decoy, decoy.indexOf('S'), CH2, 5),
				3, peaks(PEPTIDE, -1, CH2, 6), 4, peaks(PEPTIDE, S, CH2, 7)));
		List<JointCandidate> candidates = JointCandidate.merge(open, Collections.emptyList(), new int[2]);
		for (JointCandidate c : candidates) {
			c.group = JointCandidate.Group.OPEN_UNEXPLAINED;
		}
		boolean[] winner = { true, true, true, false };
		JointIdentification.Config cfg = new JointIdentification.Config();
		cfg.shiftLocaliser = new ShiftLocaliser(20);
		cfg.localisationSpectra.put("run1", mzml);
		List<String> messages = new ArrayList<>();
		JointIdentification.localiseWinners(candidates, winner, cfg, messages::add, messages::add);
		// CH2 on S: S to T has that mass, so the fragments cannot decide.
		assertEquals(JointCandidate.Group.OPEN_SUBSTITUTION_SITE, candidates.get(0).getGroup());
		assertEquals(JointCandidate.Group.OPEN_SUBSTITUTION_SITE, candidates.get(1).getGroup());
		assertEquals("OPEN_SUBSTITUTION_SITE", candidates.get(0).openPrimary().getIdentificationGroup());
		assertTrue(candidates.get(0).openPrimary().getShiftLocalisation() == candidates.get(0).getLocalisation(),
				"the PSM carries the localisation for the pangenome variant check");
		assertEquals(JointCandidate.Group.OPEN_UNEXPLAINED, candidates.get(2).getGroup());
		assertEquals(ShiftLocaliser.State.UNSUPPORTED, candidates.get(2).getLocalisation().state);
		assertNull(candidates.get(3).getLocalisation());
		assertEquals(JointCandidate.Group.OPEN_UNEXPLAINED, candidates.get(3).getGroup());
		assertTrue(JointIdentification.identifies(JointCandidate.Group.OPEN_LOCALISED, cfg));
		assertFalse(JointIdentification.identifies(JointCandidate.Group.OPEN_SUBSTITUTION_SITE, cfg));
		cfg.acceptUnexplained = true;
		assertFalse(JointIdentification.identifies(JointCandidate.Group.OPEN_SUBSTITUTION_SITE, cfg),
				"accepting unexplained shifts does not accept these");
		// Its own FDR group: the two groups' q-values are walked apart.
		double[] scores = { 5, 4, 3, 0 };
		double[] q = JointIdentification.groupQValues(candidates, scores, winner, false, 1.0);
		assertEquals(Math.min(1.0, 2.0 / 1), q[1], 1e-12, "the localised group's own walk (1 decoy, 1 target)");
		assertEquals(1.0 / 1, q[2], 1e-12, "the unexplained group alone (0 decoys, 1 target)");
	}

	@Test
	@DisplayName("Ions the tolerance cannot resolve from the unshifted ion's 13C peaks are no evidence: an unmodified "
			+ "spectrum whose large fragments show only their 13C peak does not localise +0.984, a deamidated one does")
	void unresolvableIonsAreNoEvidence() {
		ShiftLocaliser loc = new ShiftLocaliser(20);
		double deam = 0.984016;
		// 19.3 mDa from the +1 isotope step: unresolvable at 20 ppm above m/z 967 (charge 1), 483 (charge 2)
		assertFalse(loc.unresolvable(500 - ShiftLocaliser.PROTON, deam, 1, 500 + deam));
		assertTrue(loc.unresolvable(1200 - ShiftLocaliser.PROTON, deam, 1, 1200 + deam));
		assertTrue(loc.unresolvable(2 * 600 - 2 * ShiftLocaliser.PROTON, deam, 2, 600 + deam / 2));
		assertTrue(loc.unresolvable(800 - ShiftLocaliser.PROTON, 1.0034, 1, 800 + 1.0034), "an isotope step itself");
		assertFalse(loc.unresolvable(1200 - ShiftLocaliser.PROTON, CH2, 1, 1200 + CH2));
		String peptide = "GLNAVEKTPLSNYEDAWKR";
		int n = peptide.indexOf('N');
		// Above m/z 1000 the peak list holds a fragment's 13C peak and not its monoisotopic one, as real
		// peak lists often do for large fragments; without the rule this reads as +0.984 on the backbone.
		ShiftLocaliser.Result unmodified = loc.localise(peptide, 2, deam, isotopeOnlyAbove(peaks(peptide, -1, 0, 11), 1000));
		assertEquals(ShiftLocaliser.State.UNSUPPORTED, unmodified.state,
				"13C peaks of an unmodified peptide are not deamidation: support " + unmodified.support);
		ShiftLocaliser.Result deamidated = loc.localise(peptide, 2, deam, withIsotopes(peaks(peptide, n, deam, 12)));
		assertTrue(deamidated.support >= ShiftLocaliser.DEFAULT_MIN_SUPPORT,
				"the resolvable ions still carry the shift: support " + deamidated.support);
	}

	/** {@code p} with every b/y peak above {@code mz} replaced by its +1 13C peak. */
	private static MzmlScanReader.Peaks isotopeOnlyAbove(MzmlScanReader.Peaks p, double mz) {
		double[] out = p.mz.clone();
		for (int i = 0; i < out.length; i++) {
			if (p.intensity[i] > 1000f && out[i] > mz) {
				out[i] += C13;
			}
		}
		Arrays.sort(out);
		float[] in = new float[out.length];
		Arrays.fill(in, 1000f);
		return new MzmlScanReader.Peaks(out, in);
	}

	private static final double C13 = 1.0033548;

	/** {@code p} plus each b/y peak's +1 13C peak at 60% intensity, as a real spectrum has. */
	private static MzmlScanReader.Peaks withIsotopes(MzmlScanReader.Peaks p) {
		TreeMap<Double, Float> spec = new TreeMap<>();
		for (int i = 0; i < p.mz.length; i++) {
			spec.put(p.mz[i], p.intensity[i]);
			if (p.intensity[i] > 1000f) {
				spec.put(p.mz[i] + C13, p.intensity[i] * 0.6f);
			}
		}
		double[] mz = new double[spec.size()];
		float[] in = new float[spec.size()];
		int i = 0;
		for (Map.Entry<Double, Float> e : spec.entrySet()) {
			mz[i] = e.getKey();
			in[i++] = e.getValue();
		}
		return new MzmlScanReader.Peaks(mz, in);
	}

	@Test
	@DisplayName("Peak test verdicts: bounded chance share against 10%, and no verdict without a null or a scored winner")
	void peakTestDecision() {
		double max = PeakFragmentTest.DEFAULT_MAX_CHANCE_SHARE;
		assertEquals(0.10, max, 0.0);
		assertTrue(PeakFragmentTest.decide(1, 15.995, 100, 80, 2000, 10, max).passed, "Met-ox: 0.5% / 80%");
		assertFalse(PeakFragmentTest.decide(2, 4.05, 1000, 19, 2000, 34, max).passed, "band: 1.7% / 1.9%");
		assertFalse(PeakFragmentTest.decide(3, 4.05, 100, 0, 2000, 0, max).passed, "nothing supported anywhere");
		assertFalse(PeakFragmentTest.decide(6, 4.11, 320, 1, 2000, 0, max).passed,
				"a null of 0 in 2,000 is not 'no chance': 1 supported of 320 stays at the chance level");
		assertTrue(PeakFragmentTest.decide(7, -17.03, 6366, 5809, 2000, 79, max).passed,
				"-17.027: a high null (b/y - NH3 ions) but the fragments carry it in 91%");
		PeakFragmentTest.PeakResult noNull = PeakFragmentTest.decide(4, 4.05, 100, 2, 0, 0, max);
		assertTrue(noNull.passed);
		assertTrue(noNull.reason.startsWith("not tested"), noNull.reason);
		assertEquals(PeakFragmentTest.wilsonUpper(5, 400) / 0.4, PeakFragmentTest.decide(5, 42.01, 50, 20, 400, 5, max)
				.chanceShare(), 1e-12);
		// one-sided 95% Wilson bounds against hand values
		assertEquals(0.00135, PeakFragmentTest.wilsonUpper(0, 2000), 1e-5);
		assertEquals(0.04731, PeakFragmentTest.wilsonUpper(79, 2000), 1e-5);
	}

	@Test
	@DisplayName("One rule for every shift: a peak at the chance rate is withdrawn, a real one kept, and a closed "
			+ "isotope offset in another peak is localised, even one withdrawn as an identification group, while one in "
			+ "the zero family's isotope peak is not")
	void peakTestAndClosedIsotopeOffsets() throws IOException {
		double band = 4.05;
		double deamidation = 0.984016;
		String deam = "AGNVLESTPLR";
		int n = deam.indexOf('N');
		Map<Integer, MzmlScanReader.Peaks> spectra = new TreeMap<>();
		List<OpenSearchPsm> open = new ArrayList<>();
		List<OpenSearchPsm> closed = new ArrayList<>();
		int scan = 1;
		// family 10: CH2, every winner's fragments carry it
		for (int i = 0; i < 4; i++, scan++) {
			open.add(inPeak(openPsm(scan, PEPTIDE, true, CH2), 10, 10, 0));
			spectra.put(scan, peaks(PEPTIDE, S, CH2, 100 + scan));
		}
		// family 20: the band; one winner of ten shows the shift, as one unmodified PSM of ten does at that mass
		int firstBand = scan;
		for (int i = 0; i < 10; i++, scan++) {
			open.add(inPeak(openPsm(scan, PEPTIDE, true, band), 20, 20, 0));
			spectra.put(scan, peaks(PEPTIDE, i == 0 ? S : -1, band, 200 + scan));
		}
		// family 30: deamidation, seen by the open search (defines the peak's mass range)
		int openDeam = scan;
		for (int i = 0; i < 3; i++, scan++) {
			open.add(inPeak(openPsm(scan, deam, true, deamidation + (i - 1) * 0.002), 30, 30, 0));
			spectra.put(scan, peaks(deam, n, deamidation, 300 + scan));
		}
		// the zero family (0) and its +1 isotope peak (observed 31)
		int openZero = scan;
		for (int i = 0; i < 2; i++, scan++) {
			open.add(inPeak(openPsm(scan, deam, true, 1.0034 + (i - 0.5) * 0.004), 31, 0, 1));
			spectra.put(scan, peaks(deam, -1, 0, 400 + scan));
		}
		// closed winners: deamidated at isotope offset 1, a true isotope error, and 60 unmodified (the null),
		// one of which shows the band mass by chance
		int closedDeam = scan++;
		closed.add(closedPsm(closedDeam, deam, 0.9840));
		spectra.put(closedDeam, peaks(deam, n, deamidation, 501));
		int closedIso = scan++;
		closed.add(closedPsm(closedIso, deam, 1.0034));
		spectra.put(closedIso, peaks(deam, -1, 0, 502));
		int firstNull = scan;
		for (int i = 0; i < 60; i++, scan++) {
			closed.add(closedPsm(scan, PEPTIDE, 0.0005));
			spectra.put(scan, peaks(PEPTIDE, i == 0 ? S : -1, band, 600 + scan));
		}
		File mzml = new File(tmp, "run1.mzML");
		writeMzml(mzml, spectra);

		JointIdentification.prepareClosed(closed);
		List<JointCandidate> candidates = JointCandidate.merge(open, closed, new int[2]);
		boolean[] winner = new boolean[candidates.size()];
		for (int i = 0; i < candidates.size(); i++) {
			JointCandidate c = candidates.get(i);
			winner[i] = c.getScan() < openDeam || c.getScan() >= closedDeam;
			c.group = c.closed != null ? JointCandidate.Group.EXACT : JointCandidate.Group.OPEN_ACCEPTED_SHIFT;
		}
		List<ModificationValidator.PeakVerdict> verdicts = List.of(verdict(0, 0.0), verdict(10, CH2), verdict(20, band),
				verdict(30, deamidation));
		// the deamidation peak (30) is fitted but withdrawn, as substitution-shaped peaks are
		ModificationValidator.Gate gate = new ModificationValidator.Gate(new java.util.HashSet<>(List.of(0, 10, 20)),
				verdicts, java.util.Set.of());
		Map<Integer, ModificationValidator.MassClass> classes = Map.of(0, ModificationValidator.MassClass.ZERO, 10,
				ModificationValidator.MassClass.ESTABLISHED, 20, ModificationValidator.MassClass.ESTABLISHED, 30,
				ModificationValidator.MassClass.ESTABLISHED);
		JointIdentification.Config cfg = new JointIdentification.Config();
		cfg.shiftLocaliser = new ShiftLocaliser(20);
		cfg.localisationSpectra.put("run1", mzml);
		List<String> messages = new ArrayList<>();
		List<PeakFragmentTest.PeakResult> tests = JointIdentification.localiseWinners(candidates, winner, cfg,
				JointIdentification.PeakContext.of(open, gate, classes), messages::add, messages::add);

		Map<Integer, PeakFragmentTest.PeakResult> byPeak = new TreeMap<>();
		for (PeakFragmentTest.PeakResult r : tests) {
			byPeak.put(r.peakId, r);
		}
		assertTrue(byPeak.get(10).passed, "CH2 peak: " + byPeak.get(10).reason);
		assertFalse(byPeak.get(20).passed, "band at the chance rate: " + byPeak.get(20).reason);
		assertEquals(60, byPeak.get(20).nullScored);
		assertTrue(byPeak.get(30).passed, "deamidation peak: " + byPeak.get(30).reason);
		assertEquals(1, byPeak.get(30).scored, "only the closed winner is in the deamidation peak");

		JointCandidate ch2 = candidate(candidates, 1);
		assertEquals(ShiftLocaliser.State.LOCALISED_EXPLAINED, ch2.getLocalisation().state);
		JointCandidate bandWinner = candidate(candidates, firstBand);
		assertEquals(ShiftLocaliser.State.UNSUPPORTED, bandWinner.getLocalisation().state,
				"the one band PSM supported by chance is withdrawn with its peak");
		assertTrue(bandWinner.getLocalisation().explanations.isEmpty());
		assertEquals(JointCandidate.Group.OPEN_UNEXPLAINED, bandWinner.getGroup(),
				"UNIMOD has no modification of the band mass, so its unplaced winners are no identifications");

		JointCandidate deamWinner = candidate(candidates, closedDeam);
		assertEquals(JointCandidate.Group.EXACT, deamWinner.getGroup());
		OpenSearchPsm deamPsm = deamWinner.closed;
		assertEquals(ShiftLocaliser.State.LOCALISED_EXPLAINED, deamPsm.getShiftLocalisation().state);
		assertEquals(deamidation, deamPsm.getShiftLocalisation().explainedMass, 1e-4);
		assertTrue(deamPsm.isShiftLocalisedExplained());
		assertEquals(deam + "[+0.9840]", deamPsm.getModificationFamilyKey());
		assertEquals(deamPsm.getCalcNeutralPepMass() + deamPsm.getShiftLocalisation().explainedMass,
				deamPsm.getFlashLfqMonoisotopicMass(), 1e-9, "FlashLFQ looks for the deamidated form");

		JointCandidate isoWinner = candidate(candidates, closedIso);
		assertNull(isoWinner.getLocalisation(), "a +1 isotope of the zero family is an isotope error, not a shift");
		assertNull(isoWinner.closed.getShiftLocalisation());
		assertEquals(deam, isoWinner.closed.getModificationFamilyKey());
		assertNull(candidate(candidates, firstNull).getLocalisation(), "null PSMs are only scored, never labelled");
		assertNull(candidate(candidates, openDeam).getLocalisation(), "not a winner");
		assertNull(candidate(candidates, openZero).getLocalisation(), "not a winner");
		assertTrue(messages.stream().anyMatch(m -> m.contains("1 of 3 peaks not on the fragments")),
				messages.toString());
		assertTrue(messages.stream().anyMatch(m -> m.contains("2 exact-mass winners with an isotope offset: 1 in an "
				+ "isotope peak of the zero family, 1 in another fitted peak")), messages.toString());

		// r1 counts accepted winners: with none accepted, no peak can be judged and nothing is withdrawn
		for (JointCandidate c : candidates) {
			c.localisation = null;
			c.group = c.closed != null ? JointCandidate.Group.EXACT : JointCandidate.Group.OPEN_ACCEPTED_SHIFT;
		}
		JointIdentification.Localisation again = JointIdentification.localise(candidates, winner, cfg,
				JointIdentification.PeakContext.of(open, gate, classes), messages::add, messages::add);
		for (PeakFragmentTest.PeakResult r : again.testPeaks(new boolean[candidates.size()], messages::add)) {
			assertTrue(r.passed && r.reason.startsWith("not tested"), r.peakId + ": " + r.reason);
			assertEquals(0, r.scored);
		}
		ShiftLocaliser.State own = candidate(candidates, firstBand).getLocalisation().state;
		assertTrue(own == ShiftLocaliser.State.LOCALISED_EXPLAINED || own == ShiftLocaliser.State.LOCALISED_UNEXPLAINED,
				"untested, the band PSM keeps its own localisation: " + own);
	}

	// ---------------------------------------------------------------- helpers

	private static OpenSearchPsm inPeak(OpenSearchPsm p, int observed, int family, int isotopeOffset) {
		p.setObservedPeakId(observed);
		p.setPeakId(family);
		p.setPeakIsotopeOffset(isotopeOffset);
		p.setNormalizedMassDiff(p.getMassDiff() - isotopeOffset * JointIdentification.ISOTOPE_SPACING);
		return p;
	}

	private static OpenSearchPsm closedPsm(int scan, String peptide, double dm) {
		OpenSearchPsm p = openPsm(scan, peptide, true, dm);
		p.setClosedStage(true);
		p.setCalcNeutralPepMass(mass(peptide) + ShiftLocaliser.H2O);
		return p;
	}

	private static ModificationValidator.PeakVerdict verdict(int peakId, double centroid) {
		return new ModificationValidator.PeakVerdict(peakId, centroid,
				peakId == 0 ? ModificationValidator.MassClass.ZERO : ModificationValidator.MassClass.ESTABLISHED,
				100, 1, 0.01, 0.02, 0.1, false, true, "test", 0);
	}

	private static JointCandidate candidate(List<JointCandidate> candidates, int scan) {
		for (JointCandidate c : candidates) {
			if (c.getScan() == scan) {
				return c;
			}
		}
		throw new AssertionError("no candidate for scan " + scan);
	}

	private static double mass(String residues) {
		double m = 0;
		for (char c : residues.toCharArray()) {
			m += SubstitutionMassTable.residueMass(c);
		}
		return m;
	}

	private static OpenSearchPsm openPsm(int scan, String peptide, boolean target, double dm) {
		OpenSearchPsm p = new OpenSearchPsm();
		p.setFileName("run1");
		p.setScan(scan);
		p.setPeptide(peptide);
		p.setTarget(target);
		p.setHitRank(1);
		p.setHyperscore(20);
		p.setCharge(2);
		p.setMassDiff(dm);
		p.setNormalizedMassDiff(dm);
		p.setPeakId(-1);
		p.setClosedStage(false);
		p.setProtein((target ? "" : "rev_") + "MGYG000000001_1");
		return p;
	}

	/** b/y ions (charge 1) of {@code peptide}, the shift on residue {@code site} (-1: none), 85% present, plus noise. */
	private static MzmlScanReader.Peaks peaks(String peptide, int site, double delta, long seed) {
		Random r = new Random(seed);
		int n = peptide.length();
		double[] prefix = new double[n + 1];
		for (int i = 0; i < n; i++) {
			prefix[i + 1] = prefix[i] + SubstitutionMassTable.residueMass(peptide.charAt(i));
		}
		TreeMap<Double, Float> spec = new TreeMap<>();
		for (int j = 1; j <= n - 1; j++) {
			double b = prefix[j] + (site >= 0 && ShiftLocaliser.carriesSite(true, j, site, n) ? delta : 0);
			double y = prefix[n] - prefix[n - j] + ShiftLocaliser.H2O
					+ (site >= 0 && ShiftLocaliser.carriesSite(false, j, site, n) ? delta : 0);
			if (r.nextDouble() < 0.85) {
				spec.put(b + ShiftLocaliser.PROTON, 2000f);
			}
			if (r.nextDouble() < 0.85) {
				spec.put(y + ShiftLocaliser.PROTON, 3000f);
			}
		}
		for (int i = 0; i < 120; i++) {
			spec.put(150 + r.nextDouble() * 1300, 100f);
		}
		double[] mz = new double[spec.size()];
		float[] in = new float[spec.size()];
		int i = 0;
		for (Map.Entry<Double, Float> e : spec.entrySet()) {
			mz[i] = e.getKey();
			in[i++] = e.getValue();
		}
		return new MzmlScanReader.Peaks(mz, in);
	}

	private static void writeMzml(File f, Map<Integer, MzmlScanReader.Peaks> scans) throws IOException {
		StringBuilder xml = new StringBuilder("<?xml version=\"1.0\" encoding=\"utf-8\"?>"
				+ "<indexedmzML xmlns=\"http://psi.hupo.org/ms/mzml\"><mzML><run id=\"r\"><spectrumList count=\""
				+ scans.size() + "\">");
		int index = 0;
		for (Map.Entry<Integer, MzmlScanReader.Peaks> e : new TreeMap<>(scans).entrySet()) {
			double[] mz = e.getValue().mz;
			double[] in = new double[mz.length];
			for (int i = 0; i < in.length; i++) {
				in[i] = e.getValue().intensity[i];
			}
			xml.append("<spectrum index=\"").append(index++).append("\" id=\"scan=").append(e.getKey())
					.append("\" defaultArrayLength=\"").append(mz.length).append("\">")
					.append("<cvParam cvRef=\"MS\" accession=\"MS:1000511\" name=\"ms level\" value=\"2\"/>")
					.append("<binaryDataArrayList count=\"2\">").append(array("MS:1000514", mz))
					.append(array("MS:1000515", in)).append("</binaryDataArrayList></spectrum>");
		}
		xml.append("</spectrumList></run></mzML></indexedmzML>");
		Files.write(f.toPath(), xml.toString().getBytes(StandardCharsets.UTF_8));
	}

	private static String array(String kind, double[] values) {
		ByteBuffer b = ByteBuffer.allocate(values.length * 8).order(ByteOrder.LITTLE_ENDIAN);
		Arrays.stream(values).forEach(b::putDouble);
		return "<binaryDataArray encodedLength=\"0\"><cvParam cvRef=\"MS\" accession=\"MS:1000523\"/>"
				+ "<cvParam cvRef=\"MS\" accession=\"MS:1000576\"/><cvParam cvRef=\"MS\" accession=\"" + kind + "\"/>"
				+ "<binary>" + Base64.getEncoder().encodeToString(b.array()) + "</binary></binaryDataArray>";
	}
}
