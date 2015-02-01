package org.uma.jmetal.algorithm.totrythemeasures;

import org.uma.jmetal.algorithm.multiobjective.nsgaii.NSGAII;
import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.SelectionOperator;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.util.evaluator.SolutionListEvaluator;
import org.uma.jmetal.util.measurement.impl.CountingMeasure;

/**
 * Created by Antonio J. Nebro on 30/10/14.
 */
public class NSGAIIM extends NSGAII {
  private CountingMeasure iterations ;

  /**
   * Constructor
   */
  public NSGAIIM(Problem problem, int maxIterations, int populationSize,
      CrossoverOperator crossoverOperator, MutationOperator mutationOperator,
      SelectionOperator selectionOperator, SolutionListEvaluator evaluator) {
    super(problem, maxIterations, populationSize, crossoverOperator, mutationOperator, selectionOperator, evaluator) ;

  }

  @Override protected void initProgress() {
    iterations = new CountingMeasure(1) ;
  }

  @Override protected void updateProgress() {
    iterations.increment();
  }

  @Override protected boolean isStoppingConditionReached() {
    return iterations.get() >= maxIterations;
  }

}
