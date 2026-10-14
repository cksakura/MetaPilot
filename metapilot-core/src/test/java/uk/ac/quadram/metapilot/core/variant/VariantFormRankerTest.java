package uk.ac.quadram.metapilot.core.variant;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

class VariantFormRankerTest {

	private static final double G_TO_A = 71.03711 - 57.02146;

	/** GASDK with G1->A: the same +14.016 also fits S3->T and D4->E. */
	private static VariantCandidate candidate() {
		VariantCandidate c = new VariantCandidate();
		c.engine = VariantCandidate.Engine.OPEN_SEARCH;
		c.basePeptide = "GASDK";
		c.position = 1;
		c.fromAa = 'G';
		c.toAa = 'A';
		c.variantPeptide = "AASDK";
		c.deltaMass = G_TO_A;
		c.charge = 2;
		return c;
	}

	/** b2, b3, y1, y2, y3, all singly charged. */
	private static final int[] IONS = { 2, 3, 1, 2, 3 };
	private static final boolean[] Y = { false, false, true, true, true };

	private static DiannFragmentPredictor.Prediction prediction(float... intensity) {
		return new DiannFragmentPredictor.Prediction(IONS.clone(), Y.clone(), new int[] { 1, 1, 1, 1, 1 },
				intensity);
	}

	/** A spectrum with the given intensities at the fragments of {@code backbone} (shift on site 0). */
	private static MzmlScanReader.Peaks spectrum(String backbone, int shiftSite, double shift, float... intensity) {
		TreeMap<Double, Float> peaks = new TreeMap<>();
		for (int i = 0; i < IONS.length; i++) {
			peaks.put(VariantFormRanker.fragmentMz(backbone, IONS[i], Y[i], 1, shiftSite, shift), intensity[i]);
		}
		peaks.put(1000.0, 5f);
		double[] mz = new double[peaks.size()];
		float[] in = new float[peaks.size()];
		int i = 0;
		for (Map.Entry<Double, Float> e : peaks.entrySet()) {
			mz[i] = e.getKey();
			in[i++] = e.getValue();
		}
		return new MzmlScanReader.Peaks(mz, in);
	}

	private static VariantFormRanker.Predictions predictions(Map<String, DiannFragmentPredictor.Prediction> byBackbone) {
		return (backbone, charge) -> byBackbone.get(backbone);
	}

	@Test
	@DisplayName("b/y m/z, with the shift only on fragments that contain the site")
	void fragmentMz() {
		assertEquals(227.10263, VariantFormRanker.fragmentMz("PEPTIDEK", 2, false, 1, -1, 0), 1e-4);
		assertEquals(147.11280, VariantFormRanker.fragmentMz("PEPTIDEK", 1, true, 1, -1, 0), 1e-4);
		assertEquals((227.10263 + 14.01565 + VariantFormRanker.PROTON) / 2,
				VariantFormRanker.fragmentMz("PEPTIDEK", 2, false, 2, 0, 14.01565), 1e-4);
		assertEquals(147.11280, VariantFormRanker.fragmentMz("PEPTIDEK", 1, true, 1, 0, 14.01565), 1e-4,
				"y1 does not contain site 0");
		assertEquals(147.11280 + 14.01565, VariantFormRanker.fragmentMz("PEPTIDEK", 1, true, 1, 7, 14.01565), 1e-4);
		// C carries carbamidomethyl, as in the search and in DIA-NN's C[UniMod:4]
		assertEquals(103.00919 + 57.02146 + VariantFormRanker.PROTON, VariantFormRanker.fragmentMz("CK", 1, false, 1,
				-1, 0), 1e-4);
	}

	@Test
	@DisplayName("Entropy similarity: 1 for the same shape, 0 for disjoint or empty")
	void entropySimilarity() {
		assertEquals(1.0, VariantFormRanker.entropySimilarity(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 }, 3),
				1e-9);
		assertEquals(0.0, VariantFormRanker.entropySimilarity(new double[] { 1, 0 }, new double[] { 0, 1 }, 2), 1e-9);
		assertEquals(0.0, VariantFormRanker.entropySimilarity(new double[] { 1, 2 }, new double[] { 0, 0 }, 2));
		double partial = VariantFormRanker.entropySimilarity(new double[] { 10, 30, 80, 50, 100 },
				new double[] { 100, 50, 80, 30, 10 }, 5);
		assertTrue(partial > 0.6 && partial < 0.9, "reversed pattern: " + partial);
	}

	@Test
	@DisplayName("Forms: own substitution first, other same-mass substitutions, the delta on every residue")
	void forms() {
		List<VariantFormRanker.Form> forms = VariantFormRanker.forms(candidate(), 0.005);
		assertEquals(VariantFormRanker.Kind.SUBSTITUTION, forms.get(0).kind);
		assertEquals("AASDK", forms.get(0).backbone);
		assertEquals(0, forms.get(0).site);
		List<String> others = new ArrayList<>();
		List<Integer> otherSites = new ArrayList<>();
		int mods = 0;
		for (VariantFormRanker.Form f : forms) {
			if (f.kind == VariantFormRanker.Kind.OTHER_SUBSTITUTION) {
				others.add(f.backbone);
				otherSites.add(f.site);
			} else if (f.kind == VariantFormRanker.Kind.MODIFICATION) {
				assertEquals("GASDK", f.backbone);
				assertEquals(mods, f.shiftSite);
				assertEquals(mods, f.site);
				mods++;
			}
		}
		assertEquals(List.of("GATDK", "GASEK"), others);
		assertEquals(List.of(2, 3), otherSites);
		assertEquals(5, mods);
		String label = forms.get(forms.size() - 5).label;
		assertTrue(label.startsWith("modification +14.01") && label.endsWith("@G1"), label);

		VariantCandidate odd = candidate();
		odd.basePeptide = "GAXDK";
		assertTrue(VariantFormRanker.forms(odd, 0.005).isEmpty());
	}

	@Test
	@DisplayName("Verdicts: the spectrum's intensity pattern picks the form")
	void verdicts() {
		float[] variantPattern = { 100, 50, 80, 30, 10 };
		float[] basePattern = { 10, 30, 80, 50, 100 };
		// same m/z for the substitution and the modification at G1; only intensities differ
		MzmlScanReader.Peaks looksLikeVariant = spectrum("AASDK", -1, 0, variantPattern);
		MzmlScanReader.Peaks looksLikeModified = spectrum("AASDK", -1, 0, basePattern);

		Map<String, DiannFragmentPredictor.Prediction> p = new HashMap<>();
		p.put("AASDK", prediction(variantPattern));
		p.put("GASDK", prediction(basePattern));

		VariantCandidate c = candidate();
		VariantFormRanker.rank(c, VariantFormRanker.forms(c, 0.005), predictions(p), looksLikeVariant, 20, 0.05);
		assertEquals(VariantCandidate.Ms2Verdict.SUBSTITUTION, c.ms2Verdict);
		assertEquals(1.0, c.ms2ScoreSubstitution, 1e-6);
		assertTrue(c.ms2ScoreModification < 0.9);
		assertTrue(Double.isNaN(c.ms2ScoreOtherSubstitution), "no same-site other substitution");
		assertEquals(Boolean.FALSE, c.ms2LocalisesElsewhere);
		assertEquals("substitution G1A", c.ms2BestForm);

		c = candidate();
		VariantFormRanker.rank(c, VariantFormRanker.forms(c, 0.005), predictions(p), looksLikeModified, 20, 0.05);
		assertEquals(VariantCandidate.Ms2Verdict.MODIFICATION, c.ms2Verdict);
		assertTrue(c.ms2BestForm.startsWith("modification") && c.ms2BestForm.endsWith("@G1"), c.ms2BestForm);

		p.put("GASDK", prediction(variantPattern));
		c = candidate();
		VariantFormRanker.rank(c, VariantFormRanker.forms(c, 0.005), predictions(p), looksLikeVariant, 20, 0.05);
		assertEquals(VariantCandidate.Ms2Verdict.AMBIGUOUS, c.ms2Verdict, "equal predictions cannot separate them");

		p.remove("AASDK");
		c = candidate();
		VariantFormRanker.rank(c, VariantFormRanker.forms(c, 0.005), predictions(p), looksLikeVariant, 20, 0.05);
		assertEquals(VariantCandidate.Ms2Verdict.NOT_RANKED, c.ms2Verdict, "no prediction for the substitution");
		assertEquals(null, c.ms2LocalisesElsewhere);
	}

	@Test
	@DisplayName("A better form at another site sets the elsewhere flag, not the verdict")
	void elsewhere() {
		float[] variantPattern = { 100, 50, 80, 30, 10 };
		float[] basePattern = { 10, 30, 80, 50, 100 };
		Map<String, DiannFragmentPredictor.Prediction> p = new HashMap<>();
		p.put("AASDK", prediction(variantPattern));
		p.put("GASDK", prediction(basePattern));
		p.put("GATDK", prediction(variantPattern));

		VariantCandidate c = candidate();
		VariantFormRanker.rank(c, VariantFormRanker.forms(c, 0.005), predictions(p),
				spectrum("GATDK", -1, 0, variantPattern), 20, 0.05);
		assertEquals(1.0, c.ms2ScoreElsewhere, 1e-6);
		assertEquals(Boolean.TRUE, c.ms2LocalisesElsewhere);
		assertEquals("substitution S3T", c.ms2BestForm);
		assertTrue(c.ms2Verdict != VariantCandidate.Ms2Verdict.OTHER_SUBSTITUTION,
				"S3T is not at the candidate's site");

		c = candidate();
		VariantFormRanker.rank(c, VariantFormRanker.forms(c, 0.005), predictions(p),
				spectrum("GASDK", 2, G_TO_A, basePattern), 20, 0.05);
		assertEquals(Boolean.TRUE, c.ms2LocalisesElsewhere);
		assertTrue(c.ms2BestForm.endsWith("@S3"), c.ms2BestForm);
		assertTrue(c.ms2ScoreModification < 1.0, "the modification score is the one at G1");
	}

	@Test
	@DisplayName("Fragments above the precursor charge are not scored")
	void fragmentChargeCapped() {
		DiannFragmentPredictor.Prediction p = new DiannFragmentPredictor.Prediction(new int[] { 2 },
				new boolean[] { false }, new int[] { 3 }, new float[] { 100 });
		VariantFormRanker.Form f = new VariantFormRanker.Form(VariantFormRanker.Kind.SUBSTITUTION, "AASDK", 0, -1,
				0, "x");
		assertTrue(Double.isNaN(VariantFormRanker.score(f, p, spectrum("AASDK", -1, 0, 1, 1, 1, 1, 1), 2, 20)));
	}
}
