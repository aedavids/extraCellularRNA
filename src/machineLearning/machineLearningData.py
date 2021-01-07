'''
Created on Nov 9, 2020

@author: andrewdavidson
'''

###############################################################################
class MachineLearningData(object):
    '''
    main purpose: keep track of data meta data. I.E. source , pedigree, annotations, ..

    access data member directly
    self.name = ""
    self.source = ""
    self.notes = ""
    '''

    ###############################################################################
    def __init__(self, XTrain, yTrain,  XVal, yVal, XTest, yTest):
        '''
        data should be numpy array like 
        '''
        self.XTrain = XTrain
        self.yTrain = yTrain
        
        # use validation sets for hyper parameter tunning
        self.XVal = XVal
        self.yVal = yVal
        
        # hold out set
        self.XTest = XTest
        self.yTest = yTest
                
        self.name = ""
        self.source = ""
        self.notes = ""
        
    ###############################################################################
    def __str__(self):
        fmt = "data set name:\n{}\n\nsource:\n{}\n\nnotes:\n{}"
        ret = fmt.format(self.name, self.source, self.notes)
        return ret
    
    ###############################################################################
    def debug(self):
        print("XTrain.shape:{} yTrain.shape:{}".format(self.XTrain.shape, self.yTrain.shape))
        print("XVal.shape:{} yVal.shape:{}".format(self.XVal .shape, self.yVal.shape))        
        print("XTest.shape:{} yTest.shape:{}".format(self.XTest.shape, self.yTest.shape))
    