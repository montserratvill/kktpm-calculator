package kktpm;

import exceptions.EvaluationException;
import exceptions.MisplacedTokensException;
import exceptions.TooManyDecimalPointsException;
import org.junit.Assert;
import org.junit.Test;
import parsing.OptimizationProblem;
import parsing.XMLParser;

import java.io.File;

import static kktpm.KKTPMCalculator.*;
import static kktpm.KKTPMCalculator.getKKTPM;

public class KKTPMCalculatorTest {

    @Test
    public void testKKTPM_ZDT1() throws Throwable {
        File problemDescriptionFile = new File(
                getClass().getClassLoader().getResource("problems/zdt1.xml").toURI());
        // Read problem
        OptimizationProblem problem = XMLParser.readXML(problemDescriptionFile);
        double[] x = {
                0.005065176101100577, 0.24228979916064308, 0.638051895030216, 0.028511432927259106,
                0.6630621930951288, 0.5195760217178331, 0.8026922643429892, 0.5914570068228723,
                0.27057779733556053, 0.5922785232042694, 0.9520708133983169, 0.38334556740659276,
                0.1459470507434577, 0.6199474827869116, 0.14599438755486827, 0.9093004421866538,
                0.7417522157525359, 0.5992397765005402, 0.511315215257181, 0.012403566406870126,
                0.6919607667236274, 0.674203585728449, 0.3380080916453785, 0.6070230481955768,
                0.31699579503823117, 0.8366015515673879, 0.6278242708580304, 0.5288788299923252,
                0.2364277742638098, 0.4367130948868102};
        //double[] w = {0.0203, 0.9798};
        problem.setVector("x", x);
        double[] ideal = {-0.001, -0.001};
        double kktpm = KKTPMCalculator.getKKTPM(
                problem,
                ideal)
                .getKktpm();
        Assert.assertEquals(0.8463335902066791, kktpm, 10e-10);
    }


    @Test
    public void testKKTPM_BNH() throws Throwable {
        File problemDescriptionFile = new File(
                getClass().getClassLoader().getResource("problems/bnh.xml").toURI());
        // Read problem
        OptimizationProblem problem = XMLParser.readXML(problemDescriptionFile);
        problem.setAllVariables(new double[]{4.9999993, 3.0});
        double[] ideal = {-0.05, -0.05};
        double kktpm = KKTPMCalculator.getKKTPM(
                problem,
                ideal,
                0.001)
                .getKktpm();
        Assert.assertEquals(0.0035766871190923443, kktpm, 1e-10);
    }

    @Test
    public void testKKTPM_OSY() throws Throwable {
        File problemDescriptionFile = new File(
                getClass().getClassLoader().getResource("problems/osy.xml").toURI());
        // Read problem
        OptimizationProblem problem = XMLParser.readXML(problemDescriptionFile);
        problem.setAllVariables(new double[]{0.1, 0.0, 1, 0.5, 5, 0.5});
        problem.setConstant("pi", Math.PI);
        double[] ideal = {-300, -0.05};
        double kktpm = KKTPMCalculator.getKKTPM(
                problem,
                ideal,
                0.001)
                .getKktpm();
        Assert.assertEquals(0.9690644947096283, kktpm, 1e-10);
    }

    @Test
    public void testUsingRawData() {
        // Raw input data
        double[] x = new double[]{1};
        double[] f = new double[]{1};
        double[] z = null;
        double[] g = new double[]{-0.5};
        double[][] jacobianF = {{2}};
        double[][] jacobianG = {{-1}};
        double kktpmDirect = getDirectKKTPM(x, f, z, g, jacobianF, jacobianG, 0.0);
        double kktpmAdjusted = getAdjustedKKTPM(x, f, z, g, jacobianF, jacobianG);
        double kktpmProjected = getProjectedKKTPM(x, f, z, g, jacobianF, jacobianG);
        double kktpm = getKKTPM(x, f, z, g, jacobianF, jacobianG);
        Assert.assertEquals(0.2469135802469135, kktpmDirect, 1e-10);
        Assert.assertEquals(0.44444444444444453, kktpmAdjusted, 1e-10);
        Assert.assertEquals(0.4049382716049383, kktpmProjected, 1e-10);
        Assert.assertEquals(0.3654320987654321, kktpm, 1e-10);
    }

    @Test
    public void testUsingProblemObject() throws
            Throwable {
        // Create an empty problem object
        OptimizationProblem problem = new OptimizationProblem();
        // Define variable(s) (name and value of each variable)
        problem.setVariable("x1", 1.0);
        // Define objective function(s)
        problem.addObjective("x1^2");
        // Define constraint(s)
        problem.addConstraint("0.5-x1");
        // Define gradient(s) of objective(s)
        problem.setObjectivePartialDerivative(0, "x1", "2*x1");
        // Define gradient(s) of constraint(s)
        problem.setConstraintPartialDerivative(0, "x1", "-1");
        // Calculate KKTPM
        double kktpmDirect = KKTPMCalculator.getDirectKKTPM(problem, null, 0.0).getKktpm();
        double kktpmAdjusted = KKTPMCalculator.getAdjustedKKTPM(problem, null).getKktpm();
        double kktpmProjected = KKTPMCalculator.getProjectedKKTPM(problem, null).getKktpm();
        double kktpm = KKTPMCalculator.getKKTPM(problem, null).getKktpm();
        Assert.assertEquals(0.2469135802469135, kktpmDirect, 1e-10);
        Assert.assertEquals(0.4444444444444445, kktpmAdjusted, 1e-10);
        Assert.assertEquals(0.4049382716049383, kktpmProjected, 1e-10);
        Assert.assertEquals(0.3654320987654321, kktpm, 1e-10);
    }
}