[![Maven Central](https://img.shields.io/maven-central/v/com.github.000haitham000/kktpm-calculator.svg?label=Maven%20Central)](https://search.maven.org/search?q=g:%22com.github.000haitham000%22%20AND%20a:%22kktpm-calculator%22)

# KKTPM Calculator Version 2.1.0
Karush Kuhn Tucker Proximity Measure (KKTPM) Calculator is a Java implementation of the multiobjective optimization
metric developed by Deb et al. in their paper,
"A computationally fast convergence measure and implementation for single-, multiple-, and many-objective optimization", 
IEEE Transactions on Emerging Topics in Computational Intelligence, no. 1(4), pp.280-293, 2017.
The library calculates direct, adjusted and projected KKTPM metrics, in addition to their average (which is the default
and what is simply called KKTPM). The library supports both single and multi-objective optimization problems.

## DEPENDENCIES
1. Tx2Ex: An open source mathematical expressions parser (Apache L2 License)
2. Apache Commons Lang3 (Apache L2 License)

## FEATURES
1. Calculates direct-KKTPM, adjusted-KKTPM, projected-KKTPM and KKTPM (default). 
2. Support for both single and multiple objectives
3. Support for both unconstrained and constrained problems 
4. The ability to format input is a human readable XML format
5. Support for a wide variety of regular expressions
6. Support for using numerical gradients
7. A large number of pre-packaged test problems
8. Reporting back the number of solution evaluations consumed during the calculations
9. All the mathematical descriptive power provided by [Tx2Ex](https://github.com/000haitham000/tx2ex)

## GETTING STARTED
As a user, you can use KKTPM Calculator in a couple of ways:
1. Define your optimization problem and let KKTPM calculator do the evaluations for you.
2. Perform evaluations manually, then feed the evaluated values to KKTPM Calculator which in turn provides you with
the required metric value (i.e. KKTPM).

### Define your optimization problem
you can define your problem in either of the following ways:
1. As an XML file
2. Programmatically

#### How to put an optimization problem in XML format 
[bnh.xml](https://github.com/000haitham000/kktpm-calculator/blob/master/src/test/resources/problems/bnh.xml)
shows how BNH test problem can be represented in XML.
```
<?xml version="1.0" encoding="UTF-8"?>
<problem xmlns="http://www.coin-laboratory.com/xml/">
    <!-- Variables -->
    <variables>
        <variable>x1</variable>
        <variable>x2</variable>
    </variables>
    <!-- Objectives -->
    <objectives>
        <objective> 
            <function>4*x1^2+4*x2^2</function>
            <gradient>
                <derivative var="x1">8*x1</derivative>
                <derivative var="x2">8*x2</derivative>
            </gradient>
        </objective>
        <objective> 
            <function>(x1-5)^2+(x2-5)^2</function>
            <gradient>
                <derivative var="x1">2*(x1-5)</derivative>
                <derivative var="x2">2*(x2-5)</derivative>
            </gradient>
        </objective>
    </objectives>
    <!-- Constraints -->
    <constraints>
        <constraint> 
            <function>((x1-5)^2+x2^2-25)/25</function>
            <gradient>
                <derivative var="x1"> 2*(x1-5)/25 </derivative>
                <derivative var="x2"> 2*x2/25 </derivative>
            </gradient>
        </constraint>
        <constraint>
            <function>(-(x1-8)^2-(x2+3)^2+7.7)/7.7</function>
            <gradient>
                <derivative var="x1"> -2*(x1-8)/7.7 </derivative>
                <derivative var="x2"> -2*(x2+3)/7.7 </derivative>
            </gradient>
        </constraint>
        <constraint>
            <function> -x1 </function>
            <gradient>
                <derivative var="x1"> -1 </derivative>
                <derivative var="x2"> 0 </derivative>
            </gradient>
        </constraint>
        <constraint> 
            <function> -x2 </function>
            <gradient>
                <derivative var="x1"> 0 </derivative>
                <derivative var="x2"> -1 </derivative>
            </gradient>
        </constraint>
        <constraint> 
            <function> x1-5 </function>
            <gradient>
                <derivative var="x1"> 1 </derivative>
                <derivative var="x2"> 0 </derivative>
            </gradient>
        </constraint>
        <constraint> 
            <function> x2-3 </function>
            <gradient>
                <derivative var="x1"> 0 </derivative>
                <derivative var="x2"> 1 </derivative>
            </gradient>
        </constraint>
    </constraints>
</problem>
```

Variables, objectives and constraints are defined. Gradients of
all objectives and constraints are provided as well. It is also possible not to include gradients. In this case,
KKTPM Calculator uses numerical gradients (as shown in [wfg1.xml](https://github.com/000haitham000/kktpm-calculator/blob/master/src/test/resources/problems/wfg1.xml) below).

```
<?xml version="1.0" encoding="UTF-8"?>
<problem xmlns="http://www.coin-laboratory.com/xml/">
    <!-- Commands -->
    <commands>
        <command>A = 0.8</command>
        <command>B = 0.75</command>
        <command>C = 0.85</command>
        <command>PI = 3.14159265359</command>
    </commands>
    <!-- Variables -->
    <variables>
        <vector size="4"
                mins="0,0,0,0"
                maxs="2,4,6,8">
            z
        </vector>
    </variables>
    <!-- Objectives -->
    <objectives>
        <objective>
            <function> (1.0/14.0)*(6 * (A + (min(0.0, floor((abs(z[3] - 0.35) / abs(floor(0.35 - z[3]) + 0.35)) - B)) * A * (B - (abs(z[3] - 0.35) / abs(floor(0.35 - z[3]) + 0.35))) / B) - min(0.0, floor(C - (abs(z[3] - 0.35) / abs(floor(0.35 - z[3]) + 0.35)))) * (1.0 - A) * ((abs(z[3] - 0.35) / abs(floor(0.35 - z[3]) + 0.35)) - C) / (1.0 - C))^2+ 8 * (A + (min(0.0, floor((abs(z[4] - 0.35) / abs(floor(0.35 - z[4]) + 0.35)) - B)) * A * (B - (abs(z[4] - 0.35) / abs(floor(0.35 - z[4]) + 0.35))) / B) - min(0.0, floor(C - (abs(z[4] - 0.35) / abs(floor(0.35 - z[4]) + 0.35)))) * (1.0 - A) * ((abs(z[4] - 0.35) / abs(floor(0.35 - z[4]) + 0.35)) - C) / (1.0 - C))^0.02)+ 2 *(1 - cos((z[1]^0.02)*PI /2))*(1 - cos((z[2]^0.02)*PI /2))</function>
            <gradient>
                <derivative var="z[1]"/>
                <derivative var="z[2]"/>
                <derivative var="z[3]"/>
                <derivative var="z[4]"/>
            </gradient>
        </objective>
        <objective>
            <function>
                (1.0/14.0)*(6 * (A + (min(0.0, floor((abs(z[3] - 0.35) / abs(floor(0.35 - z[3]) + 0.35)) - B)) * A * (B - (abs(z[3] - 0.35) / abs(floor(0.35 - z[3]) + 0.35))) / B) - min(0.0, floor(C - (abs(z[3] - 0.35) / abs(floor(0.35 - z[3]) + 0.35)))) * (1.0 - A) * ((abs(z[3] - 0.35) / abs(floor(0.35 - z[3]) + 0.35)) - C) / (1.0 - C))^0.02+ 8 * (A + (min(0.0, floor((abs(z[4] - 0.35) / abs(floor(0.35 - z[4]) + 0.35)) - B)) * A * (B - (abs(z[4] - 0.35) / abs(floor(0.35 - z[4]) + 0.35))) / B) - min(0.0, floor(C - (abs(z[4] - 0.35) / abs(floor(0.35 - z[4]) + 0.35)))) * (1.0 - A) * ((abs(z[4] - 0.35) / abs(floor(0.35 - z[4]) + 0.35)) - C) / (1.0 - C))^0.02)+ 4 *(1 - cos((z[1]^0.02)*PI /2))*(1 - sin((z[2]^0.02)*PI /2))
            </function>
            <gradient>
                <derivative var="z[1]"/>
                <derivative var="z[2]"/>
                <derivative var="z[3]"/>
                <derivative var="z[4]"/>
            </gradient>
        </objective>
        <objective>
            <function>
                (1.0/14.0)*(6 * (A + (min(0.0, floor((abs(z[3] - 0.35) / abs(floor(0.35 - z[3]) + 0.35)) - B)) * A * (B - (abs(z[3] - 0.35) / abs(floor(0.35 - z[3]) + 0.35))) / B) - min(0.0, floor(C - (abs(z[3] - 0.35) / abs(floor(0.35 - z[3]) + 0.35)))) * (1.0 - A) * ((abs(z[3] - 0.35) / abs(floor(0.35 - z[3]) + 0.35)) - C) / (1.0 - C))^0.02 + 8 * (A + (min(0.0, floor((abs(z[4] - 0.35) / abs(floor(0.35 - z[4]) + 0.35)) - B)) * A * (B - (abs(z[4] - 0.35) / abs(floor(0.35 - z[4]) + 0.35))) / B) - min(0.0, floor(C - (abs(z[4] - 0.35) / abs(floor(0.35 - z[4]) + 0.35)))) * (1.0 - A) * ((abs(z[4] - 0.35) / abs(floor(0.35 - z[4]) + 0.35)) - C) / (1.0 - C))^0.02)+ 6 *(1 - (z[1]^0.02)- (cos(10* PI* (z[1]^0.02)+ PI /2)/(10* PI)))
            </function>
            <gradient>
                <derivative var="z[1]"/>
                <derivative var="z[2]"/>
                <derivative var="z[3]"/>
                <derivative var="z[4]"/>
            </gradient>
        </objective>
    </objectives>
    <!-- Constraints -->
    <constraints>
        <!-- z[1] boxing constraints -->
        <!-- 0 <= z[1] -->
        <constraint>
            <function>
                -z[1]
            </function>
            <gradient>
                <derivative var="z[1]">-1</derivative>
                <derivative var="z[2]">0</derivative>
                <derivative var="z[3]">0</derivative>
                <derivative var="z[4]">0</derivative>
            </gradient>
        </constraint>
        <!-- z[1] <= 2 -->
        <constraint>
            <function>
                z[1] - 2
            </function>
            <gradient>
                <derivative var="z[1]">1</derivative>
                <derivative var="z[2]">0</derivative>
                <derivative var="z[3]">0</derivative>
                <derivative var="z[4]">0</derivative>
            </gradient>
        </constraint>
        <!-- z[2] boxing constraints -->
        <!-- 0 <= z[2] -->
        <constraint>
            <function>
                -z[2]
            </function>
            <gradient>
                <derivative var="z[1]">0</derivative>
                <derivative var="z[2]">-1</derivative>
                <derivative var="z[3]">0</derivative>
                <derivative var="z[4]">0</derivative>
        </gradient>
        </constraint>
        <!-- z[2] <= 4 -->
        <constraint>
            <function>
                z[2] - 4
            </function>
            <gradient>
                <derivative var="z[1]">0</derivative>
                <derivative var="z[2]">1</derivative>
                <derivative var="z[3]">0</derivative>
                <derivative var="z[4]">0</derivative>
            </gradient>
        </constraint>
        <!-- z[3] boxing constraints -->
        <!-- 0 <= z[3] -->
        <constraint>
            <function>
                -z[3]
            </function>
            <gradient>
                <derivative var="z[1]">0</derivative>
                <derivative var="z[2]">0</derivative>
                <derivative var="z[3]">-1</derivative>
                <derivative var="z[4]">0</derivative>
            </gradient>
        </constraint>
        <!-- z[3] <= 6 -->
        <constraint>
            <function>
                z[3] - 6
            </function>
            <gradient>
                <derivative var="z[1]">0</derivative>
                <derivative var="z[2]">0</derivative>
                <derivative var="z[3]">1</derivative>
                <derivative var="z[4]">0</derivative>
            </gradient>
        </constraint>
        <!-- z[4] boxing constraints -->
        <!-- 0 <= z[4] -->
        <constraint>
            <function>
                -z[4]
            </function>
            <gradient>
                <derivative var="z[1]">0</derivative>
                <derivative var="z[2]">0</derivative>
                <derivative var="z[3]">0</derivative>
                <derivative var="z[4]">-1</derivative>
            </gradient>
        </constraint>
        <!-- z[4] <= 8 -->
        <constraint>
            <function>
                z[4] - 8
            </function>
            <gradient>
                <derivative var="z[1]">0</derivative>
                <derivative var="z[2]">0</derivative>
                <derivative var="z[3]">0</derivative>
                <derivative var="z[4]">1</derivative>
            </gradient>
        </constraint>
    </constraints>
</problem>
```

#### How to define an optimization problem grammatically
Several unit tests are provided to show how an optimization problem object is constructed programmatically. Here is how it is done for
a single objective problem with one variable and one constraint:
```
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
double kktpm = KKTPMCalculator.getKKTPM(problem, null).getKktpm();
```

### Evaluate manually, then use KKTPM Calculator just for KKTPM
This is the other way of using KKTPM. Assuming that you have the values of your objective function(s) evaluations,
constraints and the jacobian matrices (gradients) of both objectives and constraints, then you can feed these values
directly to KKTPM Calculator and it will just calculate the metric for you, as follows:
```
// Manually provide these values
double[] x = new double[]{1};
double[] f = new double[]{1};
double[] z = null;
double[] g = new double[]{-0.5};
double[][] jacobianF = {{2}};
double[][] jacobianG = {{-1}};
// Calculate KKTPM
double kktpm = getKKTPM(x, f, z, g, jacobianF, jacobianG);
```

## BUILT USING
Java SE Development Kit (JDK) 8 or later
(http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html)

## VERSIONING
We use [SemVer](http://semver.org/) for versioning.

## AUTHORS
Haitham Seada (http://haithamseada.com/)

## LICENSE
This project is licensed under the Apache License Version 2.0
(http://www.apache.org/licenses/LICENSE-2.0). A simple explanation of the
license in layman's terms can be found at
(http://www.apache.org/foundation/license-faq.html#WhatDoesItMEAN).