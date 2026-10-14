package uk.ac.quadram.metapilot.task.io.qpx;

/**
 * One entry of the QPX {@code additional_scores} array. Names are snake_case,
 * prefixed with the engine that produced the score (QPX's own convention, e.g.
 * {@code sage_hyperscore}) unless the score is engine-neutral ({@code q_value}).
 *
 * @see <a href="https://qpx.quantms.org/spec/scores/">QPX scores</a>
 */
public final class QpxScore {

	public final String name;
	public final double value;
	/** true = higher is better, false = lower is better, null = unknown. */
	public final Boolean higherBetter;

	public QpxScore(String name, double value, Boolean higherBetter) {
		this.name = name;
		this.value = value;
		this.higherBetter = higherBetter;
	}

	/** Adds the score to {@code out} unless its value is NaN or infinite. */
	public static void addIfFinite(java.util.List<QpxScore> out, String name, double value, Boolean higherBetter) {
		if (!Double.isNaN(value) && !Double.isInfinite(value)) {
			out.add(new QpxScore(name, value, higherBetter));
		}
	}
}
