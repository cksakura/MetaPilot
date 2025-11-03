package uk.ac.quadram.metapilot.mdb.unipept;

public class UnipProInfo {

	private String pro;
	private String name;
	private int taxId;
	private String[] ec;
	private String[] go;

	public UnipProInfo(String pro, String name, int taxId, String[] ec, String[] go) {
		this.pro = pro;
		this.name = name;
		this.taxId = taxId;
		this.ec = ec;
		this.go = go;
	}

	public String getPro() {
		return pro;
	}

	public void setPro(String pro) {
		this.pro = pro;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public int getTaxId() {
		return taxId;
	}

	public void setTaxId(int taxId) {
		this.taxId = taxId;
	}

	public String[] getEc() {
		return ec;
	}

	public void setEc(String[] ec) {
		this.ec = ec;
	}

	public String[] getGo() {
		return go;
	}

	public void setGo(String[] go) {
		this.go = go;
	}

}
