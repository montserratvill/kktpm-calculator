package parsing;

import org.junit.Assert;
import org.junit.Test;

import java.io.File;

public class XMLParserTest {

    @Test
    public void testWFG1() throws Throwable {
        File file = new File(getClass().getClassLoader().getResource("problems/wfg1.xml").toURI());
        OptimizationProblem problem = XMLParser.readXML(file);
        problem.setObjectivePartialDerivative(0, "z[1]", "789");
        System.out.println(problem.toString());
        Assert.assertEquals(789, problem.getObjectivePartialDerivative(0, "z[1]").getDerivative(), 10e-10);
        Assert.assertEquals(1, problem.getConstraintPartialDerivative(3, "z[2]").getDerivative(), 10e-10);
    }

    @Test
    public void testSingleObjective() throws Throwable {
        File file = new File(getClass().getClassLoader().getResource("problems/sample_single_objective.xml").toURI());
        OptimizationProblem problem = XMLParser.readXML(file);
        problem.setVariable("x", 3);
        Assert.assertEquals(6, problem.getObjectivePartialDerivative(0, "x").getDerivative(), 10e-10);
        Assert.assertEquals(-1, problem.getConstraintPartialDerivative(0, "x").getDerivative(), 10e-10);
    }

    @Test
    public void testBiObjective() throws Throwable {
        File file = new File(getClass().getClassLoader().getResource("problems/sample_bi_objective.xml").toURI());
        OptimizationProblem problem = XMLParser.readXML(file);
        problem.setVariable("x1", 2);
        problem.setVariable("x2", 3);
        problem.setVariable("x3", 4);
        Assert.assertEquals(24, problem.getObjectivePartialDerivative(0, "x1").getDerivative(), 10e-10);
        Assert.assertEquals(-4, problem.getObjectivePartialDerivative(0, "x2").getDerivative(), 10e-10);
        Assert.assertEquals(-3, problem.getObjectivePartialDerivative(0, "x3").getDerivative(), 10e-10);
        Assert.assertEquals(1, problem.getObjectivePartialDerivative(1, "x1").getDerivative(), 10e-10);
        Assert.assertEquals(0, problem.getObjectivePartialDerivative(1, "x2").getDerivative(), 10e-10);
        Assert.assertEquals(16, problem.getObjectivePartialDerivative(1, "x3").getDerivative(), 10e-10);
        Assert.assertEquals(1, problem.getConstraintPartialDerivative(0, "x1").getDerivative(), 10e-10);
        Assert.assertEquals(2, problem.getConstraintPartialDerivative(0, "x2").getDerivative(), 10e-10);
        Assert.assertEquals(0, problem.getConstraintPartialDerivative(0, "x3").getDerivative(), 10e-10);
        Assert.assertEquals(40, problem.getConstraintPartialDerivative(1, "x3").getDerivative(), 10e-10);
    }

    @Test
    public void testNumericalPartialDerivative() throws Throwable {
        File file = new File(getClass().getClassLoader().getResource("problems/sample_bi_objective.xml").toURI());
        OptimizationProblem problem = XMLParser.readXML(file);
        problem.setVariable("x1", 2);
        problem.setVariable("x2", 3);
        problem.setVariable("x3", 4);
        Assert.assertEquals(
                (Math.pow(2 + problem.getDelta() + 8, 2) - Math.pow(2 + 8, 2)) / problem.getDelta(),
                problem.getConstraintPartialDerivative(1, "x1").getDerivative(),
                10e-10);
        Assert.assertEquals(
                0,
                problem.getConstraintPartialDerivative(1, "x2").getDerivative(),
                10e-10);
    }
}