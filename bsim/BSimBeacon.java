package BSimBeacon2;

import java.awt.Color;
import java.awt.Graphics2D;
import java.util.Vector;
import java.util.function.ToDoubleFunction;

import javax.vecmath.Vector3d;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.BSimTicker;
import bsim.BSimUtils;
import bsim.draw.BSimDrawer;
import bsim.draw.BSimP3DDrawer;
import bsim.export.BSimLogger;
import bsim.export.BSimPngExporter;
import bsim.particle.BSimBacterium;
import processing.core.PConstants;
import processing.core.PGraphics3D;

/**
 * {@link BeaconBacterium} bacteria are set up for chemotaxis upon sufficient uptake of lanthanides.
 * The bacteria will die at high stress levels.
 */
public class BSimBeacon {

	@Parameter(names = "-export", description = "Enable export mode")
	public static boolean export = true;

	@Parameter(names = "-export_path", description = "Directory to export to")
	public static String exportPath = "./beacon-results/";

	@Parameter(names = "-bound", description = "Simulation bound size (µm)")
	public static int boundSize = 2400;

	@Parameter(names = "-time", description = "Time to run simulation for (seconds)")
	public static int simTime = boundSize * 5;

	@Parameter(names = "-c_x", description = "X-position of chemoattractant")
	public static int cX = boundSize-1;

	@Parameter(names = "-c_y", description = "Y-position of chemoattractant")
	public static int cY = 0;

	@Parameter(names = "-c_concentration", description = "Chemoattractant concentration")
	public static double cC = 1e12;

	@Parameter(names = "-c_diffusivity", description = "Chemoattractant diffusivity (µm²/s)")
	public static double diffusivity = 1500; // µm²/s

	@Parameter(names = "-c_decay", description = "Chemoattractant decay rate")
	public static double decayRate = 0.05;

	@Parameter(names = "-population", description = "Initial number of bacteria")
	public static int population = 300;

	@Parameter(names = "-viscosity", description = "Viscosity (Pa s)")
	public static double visc = 2.7e-3; // Pa s



	/** Calculates the mean of a value related to the bacteria */
	private static double calculateMean(ToDoubleFunction<BeaconBacterium> bacAttribute, 
			Vector<BSimBacterium> bacteria) {
		double sum = 0.0;
		for (BSimBacterium b : bacteria) {
			sum += bacAttribute.applyAsDouble((BeaconBacterium) b);
		}
		return sum / bacteria.size();
	}
	
	public static void main(String[] args) {
		BSimBeacon a = new BSimBeacon();
		new JCommander(a, args);
		a.run();
	}

	/*********************************************************
	 * Simulation Definition
	 *********************************************************/
	public void run() {
		/*********************************************************
		 * Create a new simulation object and set up simulation settings
		 */
		BSim sim = new BSim();
		sim.setDt(0.01);
		sim.setTimeFormat("0.00");
		sim.setSolid(true,true,true);
		sim.setBound(boundSize,boundSize,1);
		sim.setSimulationTime(simTime);
		sim.setVisc(visc);

		/*********************************************************
		 * Set up fields
		 */
		final double c = 12e8; // molecules (used for rendering)

		// Chemoattractant field
		final BSimChemicalField cField = new BSimChemicalField(sim,
				new int[]{80,1,1}, diffusivity, decayRate);

		// Position to place chemoattractant
		final Vector3d chemoattractantPosition = new Vector3d(cX, cY, 0);

		// Lanthanide field
		final BSimChemicalField lnField = new BSimChemicalField(sim, 
				new int[] {80,80,1}, 0.0, 0.0);

		/*********************************************************
		 * Set up the bacteria
		 */
		final Vector<BSimBacterium> bacteria = new Vector<>();
		while(bacteria.size() < population) {
			// Randomly position bacteria such that they are evenly spread
			BeaconBacterium p = new BeaconBacterium(sim,
					new Vector3d(Math.random()*sim.getBound().x,
							Math.random()*sim.getBound().y,
							Math.random()*sim.getBound().z),
					lnField);
			// Chemotaxis according to chemical field strength
			p.setGoal(cField);

			if(!p.intersection(bacteria)) {
				bacteria.add(p);
			}
		}

		/*********************************************************
		 * Set the ticker, define what happens each time step
		 */
		sim.setTicker(new BSimTicker() {
			@Override
			public void tick() {
				// Calculate average positions and [Ln]
				double sumX = 0;
				double sumL = 0;
				for(BSimBacterium b : bacteria) {
					b.action();
					b.updatePosition();
					sumX += b.getPosition().getX();
					sumL += ((BeaconBacterium)b).getLnIn();
				}
				// Log
				System.out.printf("Avg X: %8.0f | Avg L: %8.0f \n", 
						sumX / bacteria.size(),
						sumL / bacteria.size());

				// Maintain chemoattractant field
				cField.addQuantity(chemoattractantPosition, cC);
				cField.update();
			}
		});

		/*********************************************************
		 * Set the drawer for the simulation
		 */
		final BSimDrawer drawer = new BSimP3DDrawer(sim, 800,600) {
			@Override
			public void boundaries() {
				p3d.noFill();
				p3d.stroke(128, 128, 255);
				p3d.pushMatrix();
				p3d.translate((float)boundCentre.x,(float)boundCentre.y,(float)boundCentre.z);
				p3d.box((float)bound.x, (float)bound.y, (float)bound.z);
				p3d.popMatrix();
				p3d.noStroke();
			}

			@Override
			public void draw(Graphics2D g) {
				p3d.beginDraw();

				if(!cameraIsInitialised){
					p3d.camera((float)bound.x*0.5f, (float)bound.y*0.5f, boundSize,
							(float)bound.x*0.5f, (float)bound.y*0.5f, 0,
							0,1,0);
					cameraIsInitialised = true;
				}

				p3d.textFont(font);
				p3d.textMode(PConstants.SCREEN);

				p3d.sphereDetail(10);
				p3d.noStroke();
				p3d.background(255, 255,255);

				scene(p3d);
				boundaries();
				time();

				p3d.endDraw();
				g.drawImage(p3d.image, 0,0, null);
			}

			@Override
			public void time() {
				p3d.fill(0);
				p3d.text(sim.getFormattedTime(), 50, 50);
			}

			@Override
			public void scene(PGraphics3D p3d) {
				p3d.ambientLight(128, 128, 128);
				p3d.directionalLight(128, 128, 128, 1, 1, -1);
				draw(cField, Color.BLUE, (float)(255/c));
				for(BSimBacterium b : bacteria) {
					draw(b, Color.RED);
				}
			}
		};
		sim.setDrawer(drawer);

		// Run simulation, export to ./beacon-results/x.csv in export mode, otherwise preview mode
		if (export) {
			String filePath = BSimUtils.generateDirectoryPath(exportPath);

			BSimLogger logger = new BSimLogger(sim, filePath + "beacon.csv") {
				// Column titles
				@Override
				public void before() {
					super.before();
					this.write("Time (seconds),"
							+ "Mean bacterium position (microns),"
							+ "Mean [Ln3+] inside bacterium");
				}
				// Data
				@Override
				public void during() {
					this.write(sim.getFormattedTime()
							+ "," + calculateMean((b) -> b.getPosition().getX(), bacteria) 
							+ "," + calculateMean((b) -> b.getLnIn(), bacteria));
				}
			};
			logger.setDt(30); // Set export time step
			sim.addExporter(logger);

			BSimPngExporter imageExporter = new BSimPngExporter(sim, drawer, filePath);
			imageExporter.setDt(30);
			sim.addExporter(imageExporter);

			sim.export();
		} else {
			sim.preview();
		}
	}
}
