'''
Created on Nov 9, 2020

@author: aedavids@ucsc.edu
'''
from keras.layers import Dense
from keras.models import Sequential
from keras.optimizers import SGD
from tensorflow.keras import regularizers

###############################################################################
class BioMarkerModel(object):
    '''
    classdocs
    '''

    ###############################################################################
    def __init__(self, mlData, alpha=0.01, lambdaPenalty=0.0, name='no-name'):
        ''' 
        todo
        '''
        self.mlData = mlData
        colIdx = 1
        self.numFeatures =self.mlData.XTrain.shape[colIdx]
        
        self.alpha = alpha # learning rate
        self.lambdaPenalty = lambdaPenalty # regularlization penality
        self.name = name
        
        self.model = None
        self.history = None
        self.testResults = None
        
    ###############################################################################
    def __str__(self):
        fmt = "model name:{} alpha:{} lambda:{}\ndata set:\n{}"
        ret = fmt.format(self.name, self.alpha, self.lambdaPenalty, self.mlData)
        return ret
        
    ###############################################################################
    def run(self):
        '''
        build model
        compiles
        fits(train and validation data sets)
        evaluates(test data set)
        '''
        
        self._buildModel()
        self._compile()
        self._fit()
        self._evaluateTestData()
        
    ###############################################################################
    def _buildModel(self):
        '''
        derived class can over ride this to implement search across different 
        parameters spaces with out having to hack base class data members
        '''
        self.model = Sequential()
        
        #https://keras.io/api/layers/core_layers/dense/
        self.model.add( Dense(1, 
                         input_shape=(self.numFeatures,), 
                         activation='sigmoid' 
                         # https://keras.io/api/layers/initializers/
                         # default is "glorot_uniform"
                         , kernel_initializer='he_uniform'
                         , kernel_regularizer=regularizers.l1(self.lambdaPenalty)
                         #, bias_regularizer=regularizers.l1(lambdaPenalty)
                         #, activity_regularizer=regularizers.l1(lambdaPenalty)
                         , name="logisticRegressionNeuron"
                         ))
        
    ###############################################################################
    def _compile(self):
        '''
        todo
        '''
        # Gradient descent (with momentum) optimizer.
        # lr is learning rate == alpha
        # momentum == Beta is a value between 0 and 1. recommned value is 0.9
        opt = SGD(lr=self.alpha, momentum=0.9)
        # metric 'accuracy'
        self.model.compile(loss='binary_crossentropy', optimizer=opt, metrics=['binary_accuracy'])
        
    ###############################################################################
    def _fit(self):
        '''
        todo
        '''
        # https://keras.io/api/models/model_training_apis/
        self.history = self.model.fit(
                            self.mlData.XTrain, self.mlData.yTrain
                            #,batch_size=60
                            ,epochs=100 
                            ,verbose=0 # Verbosity mode. 0 = silent, 1 = progress bar, 2 = one line per epoch.
                            ,validation_data=(self.mlData.XVal, self.mlData.yVal)
                             )
        
        
    ###############################################################################
    def _evaluateTestData(self):
        # evalate model on test set
        # https://www.tensorflow.org/api_docs/python/tf/keras/Model#evaluate
        # returns Scalar test loss (if the model has a single output and no metrics) or 
        # list of scalars (if the model has multiple outputs and/or metrics). The attribute 
        # model.metrics_names will give you the display labels for the scalar outputs.
        self.testResults = self.model.evaluate(self.mlData.XTest, self.mlData.yTest, verbose=0)
