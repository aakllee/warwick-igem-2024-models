package BSimBeacon2;

import javax.vecmath.Vector3d;

import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.ode.BSimOdeSolver;
import bsim.ode.BSimOdeSystem;
import bsim.particle.BSimBacterium;

public class BeaconBacterium extends BSimBacterium {
	protected BeaconGRN grn;

	protected double laThreshold = 0.5; // [La] in cell required to activate chemotaxis 
	protected double laIn;              // [La] in cell
	
	protected double[] y;                // ODE variables
	
	protected BSimChemicalField laField; // Field representing [La] in the environment

	public BeaconBacterium(BSim sim, Vector3d position, BSimChemicalField LaField) {
		super(sim, position);
		this.laIn = 1.0;
		
		this.grn = new BeaconGRN();
		this.grn.setICs(0.0, 0.0);
		this.y = this.grn.getICs();
		
	}
	
	@Override 
	public void action() {
		super.action();
		
		// Solve ODEs
		y = BSimOdeSolver.rungeKutta45(grn, sim.getTime(), y, sim.getDt());
		
		// Check stress level
	}
	
	@Override
	public double pEndRun() {
		return (goal != null && movingUpGradient()
				// Only activate chemotaxis when La has passed threshold
				&& this.getLaIn() > this.getLaThreshold()) ? pEndRunUp : pEndRunElse; 
	}
	
	public void setLaIn(double La_in) { 
		this.laIn = La_in;
	}
	public void setLaThreshold(double La_threshold) {
		this.laThreshold = La_threshold; 
	}
	public double getLaIn() { return this.laIn; }
	public double getLaThreshold() { return this.laThreshold; }
	public BeaconGRN getGRN() { return this.grn; }
	
	// ODEs representing the GRN for LutH/LanM, and their effect on stress
	class BeaconGRN implements BSimOdeSystem {
		int numEq = 2; // Number of equations
		double[] ics;  // Initial conditions
		
		// Variables 
		double m, l;
		
		// Parameters
		double a, b;
		
		@Override
		public double[] derivativeSystem(double t, double[] y) {
			double[] dy = new double[numEq];
			
			return dy;
		}

		@Override
		public int getNumEq() { return this.numEq; }

		public void setICs(double a, double b) {
			this.ics = new double[]{a, b};
		}
		@Override
		public double[] getICs() { return this.ics; }
		
	}
	
}
