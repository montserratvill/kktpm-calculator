package parsing;

import exceptions.MisplacedTokensException;
import exceptions.TooManyDecimalPointsException;
import org.junit.Assert;
import org.junit.Test;

import java.util.Iterator;
import java.util.Map;

public class OptimizationProblemTest {

    @Test(expected = IllegalArgumentException.class)
    public void testMissingVariable() throws Throwable {
        OptimizationProblem problem = getOptimizationProblem();
        problem.getVariable("x3");
    }

    @Test(expected = IllegalArgumentException.class)
    public void testPartialDerivativeOfMissingVariable() throws Throwable {
        OptimizationProblem problem = getOptimizationProblem();
        problem.getObjectivePartialDerivative(0, "x3");
    }

    @Test
    public void testOptimizationProblemCreation() throws Throwable {
        OptimizationProblem problem = getOptimizationProblem();
        // Retrieve variables
        Assert.assertEquals(10, problem.getVariable("x1"), 10e-10);
        Assert.assertEquals(2, problem.getVariable("x2"), 10e-1);
        Assert.assertEquals(3, problem.getVariable("y1"), 10e-10);
        Assert.assertEquals(7, problem.getVariable("y3"), 10e-10);
        // Retrieve objectives
        Assert.assertEquals(-20, problem.getObjective(0), 10e-10);
        Assert.assertEquals(25, problem.getObjective(1), 10e-10);
        // Retrieve constraints
        Assert.assertEquals(1, problem.getConstraint(0), 10e-10);
        Assert.assertEquals(13, problem.getConstraint(1), 10e-10);
        // Retrive partial derivatives
        // First objective
        Assert.assertEquals(2, problem.getObjectivePartialDerivative(0, "x1").getDerivative(), 10e-10);
        Assert.assertEquals(0, problem.getObjectivePartialDerivative(0, "x2").getDerivative(), 10e-10);
        Assert.assertEquals(3, problem.getObjectivePartialDerivative(0, "y1").getDerivative(), 10e-10);
        Assert.assertEquals(-14, problem.getObjectivePartialDerivative(0, "y3").getDerivative(), 10e-10);
        // Second objective
        Assert.assertEquals(2, problem.getObjectivePartialDerivative(1, "x1").getDerivative(), 10e-10);
        Assert.assertEquals(1, problem.getObjectivePartialDerivative(1, "x2").getDerivative(), 10e-10);
        Assert.assertEquals(1, problem.getObjectivePartialDerivative(1, "y1").getDerivative(), 10e-10);
        Assert.assertEquals(0, problem.getObjectivePartialDerivative(1, "y3").getDerivative(), 10e-10);
        // First constraint
        Assert.assertEquals(0, problem.getConstraintPartialDerivative(0, "x1").getDerivative(), 10e-10);
        Assert.assertEquals(12, problem.getConstraintPartialDerivative(0, "x2").getDerivative(), 10e-10);
        Assert.assertEquals(0, problem.getConstraintPartialDerivative(0, "y1").getDerivative(), 10e-10);
        Assert.assertEquals(-1, problem.getConstraintPartialDerivative(0, "y3").getDerivative(), 10e-10);
        // Second constraint
        Assert.assertEquals(1, problem.getConstraintPartialDerivative(1, "x1").getDerivative(), 10e-10);
        Assert.assertEquals(0, problem.getConstraintPartialDerivative(1, "x2").getDerivative(), 10e-10);
        Assert.assertEquals(1, problem.getConstraintPartialDerivative(1, "y1").getDerivative(), 10e-10);
        Assert.assertEquals(0, problem.getConstraintPartialDerivative(1, "y3").getDerivative(), 10e-10);
        // Change all variables at once
        problem.setAllVariables(new double[]{-1, -2, -3, -4});
        // Retrieve variables
        Assert.assertEquals(-1, problem.getVariable("x1"), 10e-10);
        Assert.assertEquals(-2, problem.getVariable("x2"), 10e-10);
        Assert.assertEquals(-3, problem.getVariable("y1"), 10e-10);
        Assert.assertEquals(-4, problem.getVariable("y3"), 10e-10);
        // Test after retreival
        Assert.assertEquals(-27, problem.getObjective(0), 10e-10);
        Assert.assertEquals(-7, problem.getObjective(1), 10e-10);
        // Change a variable (by name)
        problem.setVariable("y1", 100);
        // Add a variable (by name)
        problem.setVariable("new_variable", 77);
        // Retreive variables (using iterator)
        Iterator<Map.Entry<String, Double>> it = problem.getVariablesIterator();
        while (it.hasNext()) {
            Map.Entry<String, Double> entry = it.next();
            switch (entry.getKey()) {
                case "x1":
                    Assert.assertEquals(-1, entry.getValue(), 10e-10);
                    break;
                case "x2":
                    Assert.assertEquals(-2, entry.getValue(), 10e-10);
                    break;
                case "y1":
                    Assert.assertEquals(100, entry.getValue(), 10e-10);
                    break;
                case "y3":
                    Assert.assertEquals(-4, entry.getValue(), 10e-10);
                    break;
                case "new_variable":
                    Assert.assertEquals(77, entry.getValue(), 10e-10);
                    break;
            }
        }
        // Test after retrieval
        Assert.assertEquals(282, problem.getObjective(0), 10e-10);
        Assert.assertEquals(96, problem.getObjective(1), 10e-10);
    }

    private OptimizationProblem getOptimizationProblem() throws Throwable {
        // Create an optimization problem object
        OptimizationProblem problem = new OptimizationProblem();
        // Variables
        problem.setVariable("x1", 10);
        problem.setVariable("x2", 2);
        problem.setVariable("y1", 3);
        problem.setVariable("y3", 7);
        // Objectives
        problem.addObjective("2*x1+3*y1-y3^2");
        problem.addObjective("2*x1+x2+y1");
        // Constraints
        problem.addConstraint("x2^3-y3");
        problem.addConstraint("x1+y1");
        // Partial Derivatives (objective 1)
        problem.setObjectivePartialDerivative(0, "x1", "2");
        problem.setObjectivePartialDerivative(0, "x2", "0");
        problem.setObjectivePartialDerivative(0, "y1", "3");
        problem.setObjectivePartialDerivative(0, "y3", "0-2*y3");
        // Partial Derivatives (objective 2)
        problem.setObjectivePartialDerivative(1, "x1", "2");
        problem.setObjectivePartialDerivative(1, "x2", "1");
        problem.setObjectivePartialDerivative(1, "y1", "1");
        problem.setObjectivePartialDerivative(1, "y3", "0");
        // Partial Derivatives (constraint 1)
        problem.setConstraintPartialDerivative(0, "x1", "0");
        problem.setConstraintPartialDerivative(0, "x2", "3*x2^2");
        problem.setConstraintPartialDerivative(0, "y1", "0");
        problem.setConstraintPartialDerivative(0, "y3", "0-1");
        // Partial Derivatives (constraint 2)
        problem.setConstraintPartialDerivative(1, "x1", "1");
        problem.setConstraintPartialDerivative(1, "x2", "0");
        problem.setConstraintPartialDerivative(1, "y1", "1");
        problem.setConstraintPartialDerivative(1, "y3", "0");
        // Return the created optimization problem
        return problem;
    }

    @Test
    public void testCreateProblemProgrammatically() throws Throwable {
        // Create an optimization problem object
        OptimizationProblem problem = new OptimizationProblem();
        // Variables
        problem.setVariable("x1", 10);
        problem.setVariable("x2", 2);
        problem.setVariable("y1", 3);
        problem.setVariable("y3", 7);
        // Objectives
        problem.addObjective("2*x1+3*y1-y3^2");
        problem.addObjective("2*x1+x2+y1");
        // Constraints
        problem.addConstraint("x2^3-y3");
        problem.addConstraint("x1+y1");
        // Partial Derivatives (objective 1)
        problem.setObjectivePartialDerivative(0, "x1", "2");
        problem.setObjectivePartialDerivative(0, "x2", "0");
        problem.setObjectivePartialDerivative(0, "y1", "3");
        problem.setObjectivePartialDerivative(0, "y3", "0-2*y3");
        // Partial Derivatives (objective 2)
        problem.setObjectivePartialDerivative(1, "x1", "2");
        problem.setObjectivePartialDerivative(1, "x2", "1");
        problem.setObjectivePartialDerivative(1, "y1", "1");
        problem.setObjectivePartialDerivative(1, "y3", "0");
        // Partial Derivatives (constraint 1)
        problem.setConstraintPartialDerivative(0, "x1", "0");
        problem.setConstraintPartialDerivative(0, "x2", "3*x2^2");
        problem.setConstraintPartialDerivative(0, "y1", "0");
        problem.setConstraintPartialDerivative(0, "y3", "0-1");
        // Partial Derivatives (constraint 2)
        problem.setConstraintPartialDerivative(1, "x1", "1");
        problem.setConstraintPartialDerivative(1, "x2", "0");
        problem.setConstraintPartialDerivative(1, "y1", "1");
        problem.setConstraintPartialDerivative(1, "y3", "0");
        // Retrieve variables
        System.out.format("x1 = %f%n", problem.getVariable("x1"));
        System.out.format("x2 = %f%n", problem.getVariable("x2"));
        try {
            System.out.format("x3 = %f%n", problem.getVariable("x3"));
        } catch (Throwable ex) {
            System.out.println(ex.toString());
        }
        System.out.format("y1 = %f%n", problem.getVariable("y1"));
        System.out.format("y3 = %f%n", problem.getVariable("y3"));
        // Retrieve objectives
        System.out.format("obj(0) = %f%n", problem.getObjective(0));
        System.out.format("obj(1) = %f%n", problem.getObjective(1));
        // Retrieve constraints
        System.out.format("con(0) = %f%n", problem.getConstraint(0));
        System.out.format("con(1) = %f%n", problem.getConstraint(1));
        // Retrive partial derivatives
        // First objective
        System.out.format("PD(obj-0,x1) = %f%n", problem.getObjectivePartialDerivative(0, "x1").getDerivative());
        System.out.format("PD(obj-0,x2) = %f%n", problem.getObjectivePartialDerivative(0, "x2").getDerivative());
        try {
            System.out.format("PD(obj-0,x3) = %f%n", problem.getObjectivePartialDerivative(0, "x3").getDerivative());
        } catch (Throwable ex) {
            System.out.println(ex.toString());
        }
        System.out.format("PD(obj-0,y1) = %f%n", problem.getObjectivePartialDerivative(0, "y1").getDerivative());
        System.out.format("PD(obj-0,y3) = %f%n", problem.getObjectivePartialDerivative(0, "y3").getDerivative());
        // Second objective
        System.out.format("PD(obj-1,x1) = %f%n", problem.getObjectivePartialDerivative(1, "x1").getDerivative());
        System.out.format("PD(obj-1,x2) = %f%n", problem.getObjectivePartialDerivative(1, "x2").getDerivative());
        System.out.format("PD(obj-1,y1) = %f%n", problem.getObjectivePartialDerivative(1, "y1").getDerivative());
        System.out.format("PD(obj-1,y3) = %f%n", problem.getObjectivePartialDerivative(1, "y3").getDerivative());
        // First constraint
        System.out.format("PD(con-0,x1) = %f%n", problem.getConstraintPartialDerivative(0, "x1").getDerivative());
        System.out.format("PD(con-0,x2) = %f%n", problem.getConstraintPartialDerivative(0, "x2").getDerivative());
        System.out.format("PD(con-0,y1) = %f%n", problem.getConstraintPartialDerivative(0, "y1").getDerivative());
        System.out.format("PD(con-0,y3) = %f%n", problem.getConstraintPartialDerivative(0, "y3").getDerivative());
        // Second constraint
        System.out.format("PD(con-1,x1) = %f%n", problem.getConstraintPartialDerivative(1, "x1").getDerivative());
        System.out.format("PD(con-1,x2) = %f%n", problem.getConstraintPartialDerivative(1, "x2").getDerivative());
        System.out.format("PD(con-1,y1) = %f%n", problem.getConstraintPartialDerivative(1, "y1").getDerivative());
        System.out.format("PD(con-1,y3) = %f%n", problem.getConstraintPartialDerivative(1, "y3").getDerivative());
        // Change all variables at once
        problem.setAllVariables(new double[]{-1, -2, -3, -4});
        // Retreive variables
        System.out.format("x1 = %f%n", problem.getVariable("x1"));
        System.out.format("x2 = %f%n", problem.getVariable("x2"));
        System.out.format("y1 = %f%n", problem.getVariable("y1"));
        System.out.format("y3 = %f%n", problem.getVariable("y3"));
        // Test after retreival
        System.out.format("obj-0(after retieval) = %f%n", problem.getObjective(0));
        System.out.format("obj-1(after retieval) = %f%n", problem.getObjective(1));
        // Change a variable (by name)
        problem.setVariable("y1", 100);
        // Add a variable (by name)
        problem.setVariable("new_variable", 77);
        // Retreive variables (using iterator)
        Iterator<Map.Entry<String, Double>> it = problem.getVariablesIterator();
        while (it.hasNext()) {
            Map.Entry<String, Double> entry = it.next();
            System.out.format("%s = %f%n", entry.getKey(), entry.getValue());
        }
        // Test after retreival
        System.out.format("obj-0(after retieval) = %f%n", problem.getObjective(0));
        System.out.format("obj-1(after retieval) = %f%n", problem.getObjective(1));
    }
}