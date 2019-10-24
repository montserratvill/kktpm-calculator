/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package parsing;

import exceptions.MisplacedTokensException;
import exceptions.TooManyDecimalPointsException;

import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.EndElement;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;
import java.io.File;
import java.io.FileInputStream;

/**
 * XMLParser provides the interface of reading an optimization problem
 * information from an XML input file into an OptimizationProblem object. The
 * XML file should be following the schema found
 * <a href="https://msu.edu/~seadahai/xml/problem.xsd">here</a>. The
 * mathematical formulas are parsed using MathExpressionParser library.
 *
 * @author Haitham
 */
public class XMLParser {

    /**
     * Parses the provided XML file into an OptimizationProblem object. The
     * exact XML syntax expected by this method is the one defined by the schema
     * found at: https://msu.edu/~seadahai/xml/problem.xsd
     *
     * @param file The XML file to be parsed
     * @return An OptimizationProblem object containing all the parsed
     * information.
     * @throws IllegalArgumentException
     * @throws TooManyDecimalPointsException
     * @throws MisplacedTokensException
     * @throws Throwable
     */
    public static OptimizationProblem readXML(File file) throws
            IllegalArgumentException,
            TooManyDecimalPointsException,
            MisplacedTokensException,
            Throwable {
        // Read from the specified XML file
        XMLInputFactory xmlInputFactory = XMLInputFactory.newInstance();
        xmlInputFactory.setProperty(XMLInputFactory.IS_COALESCING, true);
        XMLEventReader xmlEventReader = xmlInputFactory.createXMLEventReader(new FileInputStream(file));
        // Create an optimization problem object
        OptimizationProblem problem = new OptimizationProblem();
        // (status) shows if the parser is currently inside an objective tag or
        // a constraint tag or none of them.
        String status = null;
        while (xmlEventReader.hasNext()) {
            XMLEvent xmlEvent = xmlEventReader.nextEvent();
            if (xmlEvent.isStartElement()) {
                StartElement startElement = xmlEvent.asStartElement();
                if (startElement.getName().getLocalPart().equals("command")) {
                    // A new variable encountered
                    xmlEvent = xmlEventReader.nextEvent();
                    problem.executeCommand(
                            xmlEvent.asCharacters().getData().trim());
                } else if (startElement.getName().getLocalPart().equals("variable")) {
                    // A new variable encountered
                    xmlEvent = xmlEventReader.nextEvent();
                    problem.setVariable(
                            xmlEvent.asCharacters().getData().trim(), 0);
                } else if (startElement.getName().getLocalPart().equals("vector")) {
                    // A new vector encountered
                    // Get vector size
                    Attribute sizeAttr = startElement.getAttributeByName(new QName("size"));
                    int size = Integer.parseInt(sizeAttr.getValue().trim());
                    // Add the new vector to the problem
                    xmlEvent = xmlEventReader.nextEvent();
                    problem.setVector(
                            xmlEvent.asCharacters().getData().trim(), new double[size]);
                } else if (startElement.getName().getLocalPart().equals("objective")) {
                    // An objective tag encountered
                    status = "objective";
                } else if (startElement.getName().getLocalPart().equals("constraint")) {
                    // A constraint tag encountered
                    status = "constraint";
                } else if (startElement.getName().getLocalPart().equals("function")) {
                    xmlEvent = xmlEventReader.nextEvent();
                    // Sete the objective/constraint
                    if (status.equals("objective")) {
                        problem.addObjective(xmlEvent.asCharacters().getData().trim());
                    } else if (status.equals("constraint")) {
                        problem.addConstraint(xmlEvent.asCharacters().getData().trim());
                    }
                } else if (startElement.getName().getLocalPart().equals("derivative")) {
                    // Get the variable name
                    Attribute idAttr = startElement.getAttributeByName(new QName("var"));
                    String varName = idAttr.getValue().trim();
                    xmlEvent = xmlEventReader.nextEvent();
                    if (status.equals("objective")) {
                        if (xmlEvent.isEndElement()) {
                            // This means that the derivative is not provided.
                            // Set the current objective derivative with respect
                            // to the current variable to null.
                            problem.setObjectivePartialDerivative(
                                    problem.getObjectivesCount() - 1,
                                    varName,
                                    null);
                        } else {
                            // Set the current objective derivative with respect to 
                            // the current variable.
                            problem.setObjectivePartialDerivative(
                                    problem.getObjectivesCount() - 1,
                                    varName,
                                    xmlEvent.asCharacters().getData().trim());
                        }
                    } else if (status.equals("constraint")) {
                        if (xmlEvent.isEndElement()) {
                            // This means that the derivative is not provided.
                            // Set the current constraint derivative with respect
                            // to the current variable to null.
                            problem.setConstraintPartialDerivative(
                                    problem.getConstraintsCount() - 1,
                                    varName,
                                    null);
                        } else {
                            // Set the current contraint derivative with respect to 
                            // the current variable.
                            problem.setConstraintPartialDerivative(
                                    problem.getConstraintsCount() - 1,
                                    varName,
                                    xmlEvent.asCharacters().getData().trim());
                        }
                    }
                }
            }
            // Clear status if not inside either an objective or a constraint
            // tag.
            if (xmlEvent.isEndElement()) {
                EndElement endElement = xmlEvent.asEndElement();
                String closingTagName = endElement.getName().getLocalPart();
                if (closingTagName.equals("objective") || closingTagName.equals("constraint")) {
                    status = null;
                }
            }
        }
        // Return the problem object containing all the parsed information.
        return problem;
    }
}
