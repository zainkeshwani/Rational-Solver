package RationalSolver;

import java.util.ArrayList;
import java.util.Collections;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Rational {

	private double a1, b1, c1;
	private double a2, b2, c2;

	public final char INFINITY = '∞';
	public final char UNION = '∪';

	public Rational(double a1, double b1, double c1, double a2, double b2, double c2) {
		this.a1 = a1;
		this.b1 = b1;
		this.c1 = c1;
		this.a2 = a2;
		this.b2 = b2;
		this.c2 = c2;
	}

	// main methods

	public ArrayList<Double> getXIntercept() {
		ArrayList<Double> solutions = new ArrayList<>();
		if(getYIntercept() == 0) {
			return solutions;
		}

		String factored = getSimplifiedRational();

		// this will get top factors
		factored = factored.substring(0, factored.indexOf("--"));

		// get only the numbers from the factor list
		Pattern patternTop = Pattern.compile("-?\\d+(\\.\\d+)?");
		Matcher matcherTop = patternTop.matcher(factored);

		while (matcherTop.find()) {
			String numberStr = matcherTop.group();
			double number = Double.parseDouble(numberStr);
			solutions.add(-1.0 * number); // we multiply by -1 to reverse the oppositeSign()
		}

		return solutions;
	}

	public double getYIntercept() {return c1 / c2;}

	public ArrayList<Double> getHolesX() {

		// storing all of the discontinuities to return later
		ArrayList<Double> holeList = new ArrayList<Double>();

		String top = this.getFactored(true);
		String bottom = this.getFactored(false);

		ArrayList<Double> factorListTop = new ArrayList<>();
		ArrayList<Double> factorListBottom = new ArrayList<>();

		// get only the numbers from the factor list
		Pattern patternTop = Pattern.compile("-?\\d+(\\.\\d+)?");
		Matcher matcherTop = patternTop.matcher(top);

		while (matcherTop.find()) {
			String numberStr = matcherTop.group();
			double number = Double.parseDouble(numberStr);
			factorListTop.add(-1.0 * number); // we multiply by -1 to reverse the oppositeSign()
		}
		// get only the numbers from the factor list
		Pattern patternBottom = Pattern.compile("-?\\d+(\\.\\d+)?");
		Matcher matcherBottom = patternBottom.matcher(bottom);

		while (matcherBottom.find()) {
			String numberStr = matcherBottom.group();
			double number = Double.parseDouble(numberStr);
			factorListBottom.add(-1.0 * number); // we multiply by -1 to reverse the oppositeSign()
		}

		// checking to see if there is any overlap in factors by sorting and comparing

		Collections.sort(factorListTop);
		Collections.sort(factorListBottom);

		// System.out.println(factorListBottom.toString()); //for testing
		// System.out.println(factorListTop.toString()); //for testing

		// get list with least factors to avoid out of bounds exception

		for (int i = 0; i < factorListTop.size(); i++) {
			// check the two lists to find the discontinuities as holes

			if (factorListBottom.indexOf(factorListTop.get(i)) != -1 && holeList.indexOf(factorListTop.get(i)) == -1) {

				holeList.add(factorListTop.get(i));
			}

		}

		return holeList;
	}

	public ArrayList<Double> getHolesY() {
		// implement a way to get the y-value from the x values of holes
		ArrayList<Double> holesY = new ArrayList<>();
		ArrayList<Double> holesX = getHolesX();
		ArrayList<Double> topVals = new ArrayList<>();
		ArrayList<Double> bottomVals = new ArrayList<>();

		String top = getSimplifiedRational().substring(0, getSimplifiedRational().indexOf("--") - 1);
		String bottom = getSimplifiedRational().substring(getSimplifiedRational().indexOf("--\n") + 3);

		Pattern patternTop = Pattern.compile("-?\\d+(\\.\\d+)?");
		Matcher matcherTop = patternTop.matcher(top);

		while (matcherTop.find()) {
			String numberStr = matcherTop.group();
			double number = Double.parseDouble(numberStr);
			topVals.add(number); // no need to multiply by -1.0, we need opposite sign here
		}

		Pattern patternBottom = Pattern.compile("-?\\d+(\\.\\d+)?");
		Matcher matcherBottom = patternBottom.matcher(bottom);

		while (matcherBottom.find()) {
			String numberStr = matcherBottom.group();
			double number = Double.parseDouble(numberStr);
			bottomVals.add(number); // no need to multiply by -1.0, we need opposite sign here
		}

		int tempValTop = 0;
		int tempValBottom = 0;

		for (int i = 0; i < topVals.size(); i++) {
			tempValTop += (holesX.get(i) + topVals.get(i));
			tempValBottom += (holesX.get(i) + bottomVals.get(i));
		}

		holesY.add((double) tempValTop / tempValBottom);

		return holesY;
	}

	public ArrayList<Double> getVA() {
		ArrayList<Double> VAList = new ArrayList<>();

		double discriminant = Math.pow(b2, 2) - 4.0 * a2 * c2;
		if (discriminant < 0) {
			// imaginary answers
			VAList.add(null);
		}

		else {
			// real answers
			if(a2!=0) {
				VAList.add((-1.0 * b2 + Math.sqrt(discriminant)) / (2.0 * a2));
				VAList.add((-1.0 * b2 - Math.sqrt(discriminant)) / (2.0 * a2));
			}
			else {
				VAList.add(-1.0*c2);
			}
		}

		ArrayList<Double> holes = this.getHolesX();
		//System.out.println("HOLES TESTING IN HOLESX(): " + holes.toString());
		
		if(holes.size() < 1) {
			return VAList;
		}

		Collections.sort(VAList);
		Collections.sort(holes);

		for (int i = VAList.size() - 1; i >= 0; i--) {
			if (holes.indexOf(VAList.get(i)) != -1) {
				VAList.remove(i);
			}
		}

		return VAList;
	}

	public double getHA() {return a1==0 ? 0 : a1 / a2;}

	public String getRationalFactored() {
		// get the whole rational factored
		return getFactored(true) + dashesNeeded() + getFactored(false); 
	}

	public String getFactored(boolean numerator) {

		// HELPER OF getRationalFactored
		// get factors of the num. or denom.
		// since only quadratics, there are only two solutions possible

		ArrayList<Double> solutions = new ArrayList<Double>();
		if (numerator) {
			// want to factor the numerator

			double discriminant = Math.pow(b1, 2) - 4.0 * a1 * c1;
			if (discriminant < 0) {
				// imaginary answers, will output NaN
				solutions.add(0.0 / 0.0);
				solutions.add(0.0 / 0.0);
			} 
			else {
				// real answers
				solutions.add((-1.0 * b1 + Math.sqrt(discriminant)) / (2.0 * a1));
				solutions.add((-1.0 * b1 - Math.sqrt(discriminant)) / (2.0 * a1));
			}
		}
		else {
			// want to factor the denominator
			double discriminant = Math.pow(b2, 2) - 4.0 * a2 * c2;
			if (discriminant < 0) {
				// imaginary answers
				solutions.add(0.0 / 0);
				solutions.add(0.0 / 0);
			} else {
				// real answers
				solutions.add((-1.0 * b2 + Math.sqrt(discriminant)) / (2.0 * a2));
				solutions.add((-1.0 * b2 - Math.sqrt(discriminant)) / (2.0 * a2));
			}
		}
		//System.out.println(solutions.toString()); //for testing
		double solution1 = solutions.get(0);
		double solution2 = solutions.get(1);
			// System.out.println(solutions.toString()); //for testing
	
		return "(x" + oppositeSign(solution1) + Math.abs(solution1) + ")" + "(x" + oppositeSign(solution2)
				+ Math.abs(solution2) + ")";
		}

	public String getSimplifiedRational() {
		String top = getFactored(true);
		String bottom = getFactored(false);

		String topTemp = new String();
		String bottomTemp = new String();

		String topFactorFinal = new String();
		String bottomFactorFinal = new String();

		for (int i = 0; i < top.length() + 1; i++) {
			// go through factors
			if (topTemp.indexOf(')') == -1) {
				// if not yet a valid factor
				topTemp += top.charAt(i);
				bottomTemp += bottom.charAt(i);

			} else {
				// we have reached a valid factor
				if (!(topTemp.equals(bottomTemp))) {
					topFactorFinal += topTemp;
					bottomFactorFinal += bottomTemp;
					topTemp = "(";
					bottomTemp = "(";
				} else {
					topTemp = "(";
					bottomTemp = "(";
				}
			}
		}
		
		String[] topSplit = topFactorFinal.split(" ");
		String[] bottomSplit = bottomFactorFinal.split(" ");
		
		String topFINALFINAL = "", bottomFINALFINAL = "";
		
		int topIdx = 0, bottomIdx = 0;
		while(topIdx < topSplit.length && bottomIdx < bottomSplit.length) { 
			if(!topSplit[topIdx].equals(bottomSplit[bottomIdx])) {
				topFINALFINAL += topSplit[topIdx]+" ";
				bottomFINALFINAL += bottomSplit[bottomIdx]+" ";
			}
			topIdx++;
			bottomIdx++;
		}

		return topFINALFINAL + dashesNeeded() + bottomFINALFINAL;
	}

	public String getDomain() {
		// Initialized
		String domain = "(-" + INFINITY + ", ";
		// get all the disconts together and sorted
		ArrayList<Double> disconts = new ArrayList<>();
		//only want to add if there are holes
		disconts.addAll(getVA());
		
		if(getHolesX().size() > 0) {
			disconts.addAll(getHolesX());
		}
		
		Collections.sort(disconts);
		//System.out.println(disconts.toString());
		//getting NaN and Infinity, so there is a problem in getVA()

		for (int i = 0; i < disconts.size(); i++) {
			domain += disconts.get(i);
			domain += ")" + UNION + "(" + disconts.get(i) + ", ";
		}
		return domain + INFINITY + ")";
	}
	
	// helper methods

	public String dashesNeeded() {
		String dashes = new String("\n--");
		for (int i = 0; i < getFactored(true).length(); i++) {
			dashes += "-";
		}
		dashes += "\n";

		return dashes;
	}

	public String oppositeSign(double num) {
		// this is a factoring helper method, to make sure the factoring signs are
		// appropriate
		return num > 0 ? "-" : "+";
	}

	// printing out (formatting for readability)

	public String printYIntercept() {
		return "(0, " + getYIntercept() + ")";
	}

	public String printXIntercept() {
		if(getXIntercept().size() < 1) return new String("None");
		String xIntercepts = new String();
		for (double i : getXIntercept()) {
			xIntercepts += "(" + i + ", 0), ";
		}
		return xIntercepts.substring(0, xIntercepts.length() - 2);
	}

	public String printHolesXY() {
		if(getHolesX().size() < 1) return new String("None");
		String holesXY = new String();
		for (int i = 0; i < getHolesX().size(); i++) {
			holesXY += "(" + getHolesX().get(i) + ", " + getHolesY().get(i) + ")";
		}
		return holesXY;
	}
	
	public String printFactoredForm() {
		String top = "";
		String bottom = "";
		if(a1!=0 && a2!=0) {
			return getSimplifiedRational();
		}
		
		if(a1==0 && b1==0) {
			top += c1;
		}
		else if(a1==0) {
			top += "(" + b1 + "x" + oppositeSign(c1) + Math.abs(c1) + ")";
		}
		
		if(a2==0) {
			bottom += "(" + b2 + "x" + oppositeSign(c2) + Math.abs(c2) + ")";
		}
		
		return top + "\n---------\n" + bottom;
	}

	
}
