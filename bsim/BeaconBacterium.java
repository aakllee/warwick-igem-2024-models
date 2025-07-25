package BSimBeacon2;

import javax.vecmath.Vector3d;

import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.ode.BSimOdeSolver;
import bsim.ode.BSimOdeSystem;
import bsim.particle.BSimBacterium;
/**
 * Class representing our engineered M. extorquens. 
 * Upon sufficient uptake of lanthanides, the bacterium activates chemotaxis.
 * Mex has only one flagellum, and so is unlikely to follow the same run-and-tumble
 * mechanics as E. coli, however we make the assumption that Mex will follow a 
 * similar net movement, as the run-reverse-flick mechanics that Mex likely follows
 * is largely very similar, involving a period of running and a period of random
 * direction change. 
 * (Park, Y., Kim, Y. and Lim, S. (2019) ‘Locomotion of a single-flagellated bacterium’, 
 * Journal of Fluid Mechanics, 859, pp. 586–612. doi:10.1017/jfm.2018.799.)
 */
public class BeaconBacterium extends BSimBacterium {
	protected BeaconGRN grn;

	protected double lnThreshold = 0.5;  // [Ln3+] in cell required to activate chemotaxis 
	protected double lnIn;               // [Ln3+] in cell
	
	protected double[] y;                // ODE variables
	
	protected BSimChemicalField laField; // Field representing [Ln3+] in the environment

	public BeaconBacterium(BSim sim, Vector3d position, BSimChemicalField lnField) {
		super(sim, position);
		this.lnIn = 1.0;
		
		
		// Parameters
		this.forceMagnitude = 1;          // pN
		this.shortTermMemoryDuration = 1; // seconds
		this.longTermMemoryDuration = 1;  // seconds
		this.sensitivity = 1;             // (molecules/(µm)³)
		
		/*// See BSimBacterium
		this.pEndRunUp = 1/1.07;          // 1/t, t = mean time to end a run up (seconds)
		this.pEndRunElse = 1/0.86;        // 1/t, t = mean time to end a run up when not moving up a gradient (seconds)
		this.pEndTumble = 1/0.14;         // 1/t, t = mean time to end a tumble (seconds)
		*/
		
		// Setup gene regulatory network
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
				&& this.getLnIn() > this.getLnThreshold()) ? 
						this.pEndRunUp : this.pEndRunElse; 
	}
	
	public void setLaIn(double lnIn) { 
		this.lnIn = lnIn;
	}
	public void setLaThreshold(double lnThreshold) {
		this.lnThreshold = lnThreshold; 
	}
	public double getLnIn() { return this.lnIn; }
	public double getLnThreshold() { return this.lnThreshold; }
	public BeaconGRN getGRN() { return this.grn; }
	
	// ODEs representing the GRN for LutH/LanM, and their effect on stress
	protected class BeaconGRN implements BSimOdeSystem {
		int numEq = 2; // Number of equations
		double[] ics;  // Initial conditions
		
		// Variables 
		double m, l;
		
		// Parameters
		double a, b;
		
		@Override
		public double[] derivativeSystem(double t, double[] y) {
			double[] dy = new double[this.numEq];
			
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
