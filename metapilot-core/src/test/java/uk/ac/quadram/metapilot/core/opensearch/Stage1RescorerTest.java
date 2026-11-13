package uk.ac.quadram.metapilot.core.opensearch;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import weka.classifiers.AbstractClassifier;
import weka.classifiers.Classifier;
import weka.core.Instance;
import weka.core.Instances;

class Stage1RescorerTest {

	@Test
	void stageOneRanksOnTheLogitWithARealRidge() {
		SemiSupervisedPsmClassifier.Config cfg = new Stage1Rescorer.Config().classifier;
		assertTrue(cfg.rankByLogit);
		assertNotSame(SemiSupervisedPsmClassifier.LOGISTIC_REGRESSION, cfg.modelFactory);
		assertTrue(cfg.modelFactory.describe().contains("ridge " + Stage1Rescorer.RIDGE), cfg.modelFactory.describe());
	}

	@Test
	void sharedClassifierDefaultsAreUnchanged() {
		SemiSupervisedPsmClassifier.Config cfg = new SemiSupervisedPsmClassifier.Config();
		assertTrue(!cfg.rankByLogit);
		assertEquals("logistic regression (weka.classifiers.functions.Logistic, standardised inputs)",
				cfg.modelFactory.describe());
	}

	/**
	 * A model confident enough that every strong candidate's target probability rounds to exactly
	 * 1.0, with a few decoys at the weak end of that block. On probabilities the block is one tie,
	 * walked in row order, so decoys early in the file reject all of it. On the logit the decoys
	 * rank below the targets they are weaker than.
	 */
	@Test
	void saturatedProbabilitiesDoNotTieDecoysWithTargets() throws Exception {
		List<OpenSearchPsm> psms = new ArrayList<>();
		int scan = 1;
		for (int i = 0; i < 5; i++) {
			psms.add(psm(scan++, false, 40 + 0.2 * i));
		}
		for (int i = 0; i < 300; i++) {
			psms.add(psm(scan++, true, 42 + 0.1 * i));
		}
		for (int i = 0; i < 3000; i++) {
			psms.add(psm(scan++, i % 2 == 0, 3.0 * i / 3000));
		}
		String[] features = { "neg_log10_expect" };

		SemiSupervisedPsmClassifier.Config cfg = new SemiSupervisedPsmClassifier.Config();
		cfg.modelFactory = SATURATING;
		SemiSupervisedPsmClassifier.Result byProbability = SemiSupervisedPsmClassifier.classify(psms, cfg, features,
				null);
		cfg.rankByLogit = true;
		SemiSupervisedPsmClassifier.Result byLogit = SemiSupervisedPsmClassifier.classify(psms, cfg, features, null);

		assertEquals(0, byProbability.acceptedTargetsAtOnePercent);
		assertEquals(300, byLogit.acceptedTargetsAtOnePercent);
	}

	private static OpenSearchPsm psm(int scan, boolean target, double negLog10Expect) {
		OpenSearchPsm psm = new OpenSearchPsm();
		psm.setFileName("run1");
		psm.setScan(scan);
		psm.setHitRank(1);
		psm.setPeptide("PEPTIDE" + scan + "K");
		psm.setProtein((target ? "MGYG_" : "rev_MGYG_") + scan);
		psm.setTarget(target);
		psm.setExpect(Math.pow(10, -negLog10Expect));
		return psm;
	}

	/** Decoy probability exp(-x): above x of about 37 the target probability is exactly 1.0. */
	private static final SemiSupervisedPsmClassifier.ModelFactory SATURATING = new SemiSupervisedPsmClassifier.ModelFactory() {
		@Override
		public Classifier create() {
			return new AbstractClassifier() {
				private static final long serialVersionUID = 1L;

				@Override
				public void buildClassifier(Instances data) {
				}

				@Override
				public double[] distributionForInstance(Instance instance) {
					double decoy = Math.exp(-instance.value(0));
					return new double[] { 1 - decoy, decoy };
				}
			};
		}

		@Override
		public String describe() {
			return "saturating stub";
		}
	};
}
