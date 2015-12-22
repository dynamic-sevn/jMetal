package org.uma.jmetal.experiment;

/**
 * Created by ajnebro on 22/12/15.
 */
import org.uma.jmetal.algorithm.Algorithm;
import org.uma.jmetal.algorithm.multiobjective.nsgaii.NSGAIIBuilder;
import org.uma.jmetal.algorithm.multiobjective.smpso.SMPSOBuilder;
import org.uma.jmetal.algorithm.multiobjective.spea2.SPEA2Builder;
import org.uma.jmetal.operator.impl.crossover.SBXCrossover;
import org.uma.jmetal.operator.impl.mutation.PolynomialMutation;
import org.uma.jmetal.problem.DoubleProblem;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.problem.multiobjective.ebes.Ebes;
import org.uma.jmetal.problem.multiobjective.zdt.*;
import org.uma.jmetal.qualityindicator.impl.*;
import org.uma.jmetal.qualityindicator.impl.hypervolume.PISAHypervolume;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.util.JMetalException;
import org.uma.jmetal.util.archive.impl.CrowdingDistanceArchive;
import org.uma.jmetal.util.experiment.ExperimentConfiguration;
import org.uma.jmetal.util.experiment.ExperimentConfigurationBuilder;
import org.uma.jmetal.util.experiment.component.*;
import org.uma.jmetal.util.experiment.util.TaggedAlgorithm;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Example of experimental study based on solving Ebes
 *
 * This experiment assumes that the reference Pareto front are known, so the names of files containing
 * them and the directory where they are located must be specified.
 *
 * Six quality indicators are used for performance assessment.
 *
 * The steps to carry out the experiment are:
 * 1. Configure the experiment
 * 2. Execute the algorithms
 * 3. Generate the reference Pareto front
 * 4. Compute que quality indicators
 * 5. Generate Latex tables reporting means and medians
 * 6. Generate Latex tables with the result of applying the Wilcoxon Rank Sum Test
 *
 * @author Antonio J. Nebro <antonio@lcc.uma.es>
 */
public class EbesStudy {
  public static void main(String[] args) throws IOException {
    if (args.length < 2) {
      throw new JMetalException("Experiment directory required as program parameter") ;
    }

    String experimentDirectory = args[0] ;

    List<Problem<DoubleSolution>> problemList = Arrays.<Problem<DoubleSolution>>asList(new Ebes()) ;

    List<TaggedAlgorithm<List<DoubleSolution>>> algorithmList = configureAlgorithmList(problemList) ;

    List<String> referenceFrontFileNames = Arrays.asList() ;

    ExperimentConfiguration<DoubleSolution, List<DoubleSolution>> configuration =
        new ExperimentConfigurationBuilder<DoubleSolution, List<DoubleSolution>>("EbesStudy")
            .setAlgorithmList(algorithmList)
            .setProblemList(problemList)
            .setExperimentBaseDirectory(experimentDirectory)
            .setOutputParetoFrontFileName("FUN")
            .setOutputParetoSetFileName("VAR")
            .setReferenceFrontDirectory("/pareto_fronts")
            .setReferenceFrontFileNames(referenceFrontFileNames)
            .setIndicatorList(Arrays.asList(
                new Epsilon<DoubleSolution>(), new Spread<DoubleSolution>(), new GenerationalDistance<DoubleSolution>(),
                new PISAHypervolume<DoubleSolution>(),
                new InvertedGenerationalDistance<DoubleSolution>(), new InvertedGenerationalDistancePlus<DoubleSolution>()))
            .setIndependentRuns(30)
            .setNumberOfCores(8)
            .build();

    new ExecuteAlgorithms<>(configuration).run();
    new GenerateReferenceParetoFront(configuration).run();
    new ComputeQualityIndicators<>(configuration).run() ;
    new GenerateLatexTablesWithStatistics(configuration).run() ;
    new GenerateWilcoxonTestTablesWithR<>(configuration).run() ;
    new GenerateFriedmanTestTables<>(configuration).run();
    new GenerateBoxplotsWithR<>(configuration).setRows(1).setColumns(1).run() ;
  }

  /**
   * The algorithm list is composed of pairs {@link Algorithm} + {@link Problem} which form part of a
   * {@link TaggedAlgorithm}, which is a decorator for class {@link Algorithm}. The {@link TaggedAlgorithm}
   * has an optional tag component, that can be set as it is shown in this example, where four variants of a
   * same algorithm are defined.
   *
   * @param problemList
   * @return
   */
  static List<TaggedAlgorithm<List<DoubleSolution>>> configureAlgorithmList(List<Problem<DoubleSolution>> problemList) {
    List<TaggedAlgorithm<List<DoubleSolution>>> algorithms = new ArrayList<>() ;
    for (int i = 0; i < problemList.size(); i++) {
      Algorithm<List<DoubleSolution>> algorithm = new NSGAIIBuilder<>(problemList.get(i), new SBXCrossover(1.0, 20.0),
          new PolynomialMutation(1.0 / problemList.get(i).getNumberOfVariables(), 20.0))
          .setMaxEvaluations(25000)
          .setPopulationSize(100)
          .build();
      algorithms.add(new TaggedAlgorithm<List<DoubleSolution>>(algorithm, "NSGAII", problemList.get(i)));
    }

    for (Problem<DoubleSolution> problem : problemList) {
      Algorithm<List<DoubleSolution>> algorithm = new SPEA2Builder<DoubleSolution>(problem, new SBXCrossover(1.0, 10.0),
          new PolynomialMutation(1.0 / problem.getNumberOfVariables(), 20.0))
          .build();
      algorithms.add(new TaggedAlgorithm<List<DoubleSolution>>(algorithm, problem));
    }

    for (Problem<DoubleSolution> problem : problemList) {
      Algorithm<List<DoubleSolution>> algorithm = new SMPSOBuilder((DoubleProblem) problem,
          new CrowdingDistanceArchive<DoubleSolution>(100))
          .build();
      algorithms.add(new TaggedAlgorithm<List<DoubleSolution>>(algorithm, problem));
    }

    return algorithms ;
  }
}
