package uk.ac.quadram.metapilot.task.dia.par;

import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

/**
 * @deprecated Use {@link MetaParameterDiaMarker} instead.
 * This class is kept as a thin alias for backward compatibility with older code/serialized configs.
 */
@Deprecated
public class MetaParameterDiaFocus extends MetaParameterDiaMarker {

	public MetaParameterDiaFocus() {
		super();
	}

	public MetaParameterDiaFocus(MetaBaseParameter metaPar, MagDbItem magDbItem, boolean librarySearch,
			boolean selfModelSearch, String[] library, boolean useModel, String[] models) {
		super(metaPar, magDbItem, librarySearch, selfModelSearch, library, useModel, models);
	}
}

