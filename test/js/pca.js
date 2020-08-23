'use strict';

Object.defineProperty(exports, '__esModule', { value: true });

var mlMatrix = require('ml-matrix');

/**
 * Creates new PCA (Principal Component Analysis) from the dataset
 * @param {Matrix} dataset - dataset or covariance matrix.
 * @param {Object} [options]
 * @param {boolean} [options.isCovarianceMatrix=false] - true if the dataset is a covariance matrix.
 * @param {string} [options.method='SVD'] - select which method to use: SVD (default), covarianceMatrirx or NIPALS.
 * @param {number} [options.nCompNIPALS=2] - number of components to be computed with NIPALS.
 * @param {boolean} [options.center=true] - should the data be centered (subtract the mean).
 * @param {boolean} [options.scale=false] - should the data be scaled (divide by the standard deviation).
 * @param {boolean} [options.ignoreZeroVariance=false] - ignore columns with zero variance if `scale` is `true`.
 * */
class PCA {
  constructor(dataset, options = {}) {
    if (dataset === true) {
      const model = options;
      this.center = model.center;
      this.scale = model.scale;
      this.means = model.means;
      this.stdevs = model.stdevs;
      this.U = mlMatrix.Matrix.checkMatrix(model.U);
      this.S = model.S;
      this.R = model.R;
      this.excludedFeatures = model.excludedFeatures || [];
      return;
    }

    dataset = new mlMatrix.Matrix(dataset);

    const {
      isCovarianceMatrix = false,
      method = 'SVD',
      nCompNIPALS = 2,
      center = true,
      scale = false,
      ignoreZeroVariance = false,
    } = options;

    this.center = center;
    this.scale = scale;
    this.means = null;
    this.stdevs = null;
    this.excludedFeatures = [];

    if (isCovarianceMatrix) {
      // User provided a covariance matrix instead of dataset.
      this._computeFromCovarianceMatrix(dataset);
      return;
    }

    this._adjust(dataset, ignoreZeroVariance);
    switch (method) {
      case 'covarianceMatrix': {
        // User provided a dataset but wants us to compute and use the covariance matrix.
        const covarianceMatrix = new mlMatrix.MatrixTransposeView(dataset)
          .mmul(dataset)
          .div(dataset.rows - 1);
        this._computeFromCovarianceMatrix(covarianceMatrix);
        break;
      }
      case 'NIPALS': {
        this._computeWithNIPALS(dataset, nCompNIPALS);
        break;
      }
      case 'SVD': {
        const svd = new mlMatrix.SVD(dataset, {
          computeLeftSingularVectors: false,
          computeRightSingularVectors: true,
          autoTranspose: true,
        });

        this.U = svd.rightSingularVectors;

        const singularValues = svd.diagonal;
        const eigenvalues = [];
        for (const singularValue of singularValues) {
          eigenvalues.push((singularValue * singularValue) / (dataset.rows - 1));
        }
        this.S = eigenvalues;
        break;
      }
      default: {
        throw new Error(`unknown method: ${method}`);
      }
    }
  }

  /**
   * Load a PCA model from JSON
   * @param {Object} model
   * @return {PCA}
   */
  static load(model) {
    if (typeof model.name !== 'string') {
      throw new TypeError('model must have a name property');
    }
    if (model.name !== 'PCA') {
      throw new RangeError(`invalid model: ${model.name}`);
    }
    return new PCA(true, model);
  }

  /**
   * Project the dataset into the PCA space
   * @param {Matrix} dataset
   * @param {Object} options
   * @return {Matrix} dataset projected in the PCA space
   */
  predict(dataset, options = {}) {
    const { nComponents = this.U.columns } = options;
    dataset = new mlMatrix.Matrix(dataset);
    if (this.center) {
      dataset.subRowVector(this.means);
      if (this.scale) {
        for (let i of this.excludedFeatures) {
          dataset.removeColumn(i);
        }
        dataset.divRowVector(this.stdevs);
      }
    }
    var predictions = dataset.mmul(this.U);
    return predictions.subMatrix(0, predictions.rows - 1, 0, nComponents - 1);
  }

  /**
   * Calculates the inverse PCA transform
   * @param {Matrix} dataset
   * @return {Matrix} dataset projected in the PCA space
   */
  invert(dataset) {
    dataset = mlMatrix.Matrix.checkMatrix(dataset);

    var inverse = dataset.mmul(this.U.transpose());

    if (this.center) {
      if (this.scale) {
        inverse.mulRowVector(this.stdevs);
      }
      inverse.addRowVector(this.means);
    }

    return inverse;
  }


  /**
   * Returns the proportion of variance for each component
   * @return {[number]}
   */
  getExplainedVariance() {
    var sum = 0;
    for (const s of this.S) {
      sum += s;
    }
    return this.S.map((value) => value / sum);
  }

  /**
   * Returns the cumulative proportion of variance
   * @return {[number]}
   */
  getCumulativeVariance() {
    var explained = this.getExplainedVariance();
    for (var i = 1; i < explained.length; i++) {
      explained[i] += explained[i - 1];
    }
    return explained;
  }

  /**
   * Returns the Eigenvectors of the covariance matrix
   * @returns {Matrix}
   */
  getEigenvectors() {
    return this.U;
  }

  /**
   * Returns the Eigenvalues (on the diagonal)
   * @returns {[number]}
   */
  getEigenvalues() {
    return this.S;
  }

  /**
   * Returns the standard deviations of the principal components
   * @returns {[number]}
   */
  getStandardDeviations() {
    return this.S.map((x) => Math.sqrt(x));
  }

  /**
   * Returns the loadings matrix
   * @return {Matrix}
   */
  getLoadings() {
    return this.U.transpose();
  }

  /**
   * Export the current model to a JSON object
   * @return {Object} model
   */
  toJSON() {
    return {
      name: 'PCA',
      center: this.center,
      scale: this.scale,
      means: this.means,
      stdevs: this.stdevs,
      U: this.U,
      S: this.S,
      excludedFeatures: this.excludedFeatures,
    };
  }

  _adjust(dataset, ignoreZeroVariance) {
    if (this.center) {
      const mean = dataset.mean('column');
      const stdevs = this.scale
        ? dataset.standardDeviation('column', { mean })
        : null;
      this.means = mean;
      dataset.subRowVector(mean);
      if (this.scale) {
        for (let i = 0; i < stdevs.length; i++) {
          if (stdevs[i] === 0) {
            if (ignoreZeroVariance) {
              dataset.removeColumn(i);
              stdevs.splice(i, 1);
              this.excludedFeatures.push(i);
              i--;
            } else {
              throw new RangeError(
                `Cannot scale the dataset (standard deviation is zero at index ${i}`,
              );
            }
          }
        }
        this.stdevs = stdevs;
        dataset.divRowVector(stdevs);
      }
    }
  }

  _computeFromCovarianceMatrix(dataset) {
    const evd = new mlMatrix.EVD(dataset, { assumeSymmetric: true });
    this.U = evd.eigenvectorMatrix;
    this.U.flipRows();
    this.S = evd.realEigenvalues;
    this.S.reverse();
  }

  _computeWithNIPALS(dataset, nCompNIPALS) {
    this.U = new mlMatrix.Matrix(nCompNIPALS, dataset.columns);
    this.S = [];

    let x = dataset;
    for (let i = 0; i < nCompNIPALS; i++) {
      let dc = new mlMatrix.NIPALS(x);

      this.U.setRow(i, dc.w.transpose());
      this.S.push(Math.pow(dc.s.get(0, 0), 2));

      x = dc.xResidual;
    }
    this.U = this.U.transpose(); // to be compatible with API
  }
}

exports.PCA = PCA;
