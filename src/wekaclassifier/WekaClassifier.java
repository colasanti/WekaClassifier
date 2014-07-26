/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package wekaclassifier;

/**
 *
 * @author ric
 */

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import weka.attributeSelection.InfoGainAttributeEval;
import weka.attributeSelection.Ranker;
import weka.classifiers.Evaluation;
import weka.classifiers.bayes.NaiveBayes;
import weka.classifiers.functions.SMO;
import weka.classifiers.meta.FilteredClassifier;
import weka.classifiers.rules.PART;
import weka.classifiers.trees.J48;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.converters.ArffLoader.ArffReader;
import weka.core.tokenizers.WordTokenizer;
import weka.filters.Filter;
import weka.filters.MultiFilter;
import weka.filters.supervised.attribute.AttributeSelection;
import weka.filters.unsupervised.attribute.Remove;
import weka.filters.unsupervised.attribute.StringToWordVector;

/**
 *
 * @author ric
 */
public class WekaClassifier {
    /**
     * Mollicutes 1
Bacilli 37
Acidobacteria 3
Chlamydiae 6
Nitrospira 1
delta/epsilon subdivisions 32
Chlorobi 2
Planctomycetia 1
Prochlorales 1
Aquificae 6
Thermoplasmata 1
Gloeobacteria 1
Bacteroidetes 17
Chloroflexi 3
Negativicutes 2
Thermotogae 1
Fusobacteriia 2
Erysipelotrichia 1
Nostocales 2
Thermoprotei 3
Gammaproteobacteria 95
Deferribacteres 1
Spirochaetia 1
Clostridia 19
ssRNA positive-strand viruses 1
Verrucomicrobia 1
Betaproteobacteria 30
Actinobacteria 37
Alphaproteobacteria 70
Deinococci 5
Oscillatoriophycideae 6
     * 
     * 
     */
    
    /**
     * Object that stores training data.
     */
    Instances trainData, testData;
    FilteredClassifier classifier;
    HashMap<String,ArrayList> hm = new HashMap<String, ArrayList>();
    HashMap<String,ArrayList> excluded = new HashMap<String, ArrayList>();
    
    private Instances initInstances(String[] classes) {
        Instances newIntstances = null;
        FastVector attributes = new FastVector(2);
        attributes.addElement(new Attribute("functions", (FastVector) null));
        FastVector classValues = new FastVector(classes.length);
        for (String e : classes) {
            classValues.addElement(e);
        }
        attributes.addElement(new Attribute("Class", classValues));
        newIntstances = new Instances("genomeFunctions", attributes, 100);
        newIntstances.setClassIndex(newIntstances.numAttributes() - 1);
        return newIntstances;
    }

    private void loadRawData(String inFileName, String myType) {
        HashMap<String, Integer> map = new HashMap<String, Integer>();
        hm.put("faculatative", new ArrayList<String>());
        hm.put("anaerobe", new ArrayList<String>());
        hm.put("aerobe", new ArrayList<String>());
        excluded.put("faculatative", new ArrayList<String>());
        excluded.put("anaerobe", new ArrayList<String>());
        excluded.put("aerobe", new ArrayList<String>());
        Pattern pattern1 = Pattern.compile("[\\s[@/|]\\s][;\\s]");
        Pattern pattern2 = Pattern.compile(";\\s");
        try {
            File inFile = new File(inFileName);
            FileReader fileReader = new FileReader(inFile);
            BufferedReader reader = new BufferedReader(fileReader);

            String line;
            while ((line = reader.readLine()) != null) {
                String delims = "\t";
                String[] tokens = line.split(delims);
                String classType = "empty";
                if (map.get(tokens[2]) == null) {
                    map.put(tokens[2], 1);
                } else {
                    Integer count = map.get(tokens[2]);
                    map.put(tokens[2], count + 1);
                }
                if ((tokens[4].equals("faculatative"))
                        || (tokens[4].equals("facultative"))) {
                    classType = "faculatative";
                }
                if ((tokens[4].equals("anaerobic"))
                        || (tokens[4].equals("anaerobe"))) {
                    classType = "anaerobe";
                }
                if ((tokens[4].equals("aerobic"))
                        || (tokens[4].equals("aerobe"))) {
                    classType = "aerobe";
                }
                /*	if ((tokens[3].equals("N"))) {
                 classType = "N";
                 }
                 if ((tokens[3].equals("P"))) {
                 classType = "P";
                 } */
                if (classType != "empty") {
                    String elements = "";
                    for (int i = 5; i < tokens.length; i++) {
                        String aHash = tokens[i];
                        Matcher matcher = pattern1.matcher(aHash);
                        if (matcher.find()) {
                            aHash = matcher.replaceAll("\t");
                        }
                        matcher = pattern2.matcher(aHash);
                        if (matcher.find()) {
                            aHash = matcher.replaceAll("\t");
                        }
                        if (myType.equals("functions")) {
                            String[] tkns = aHash.split(delims);
                            for (String s : tkns) {
                                s = s.trim();
                                s = s.replaceAll("[ '.,;]", "_");
                                if (!(s.toLowerCase().contains("ypothetical"))) {
                                    elements = elements + "\t" + s;
                                }
                            }
                        } else {
                            elements = elements + "\tH" + Integer.toString(tokens[i].hashCode());
                        }
                    }
                    //System.out.println(tokens[2]);
                    if(tokens[2].contains("Alphaproteobacteria")){
                        ArrayList a = (ArrayList) excluded.get(classType);
                       // System.out.println(tokens[2]+"*****");
                        a.add(elements);
                    }else{
                        ArrayList a = (ArrayList) hm.get(classType);
                        a.add(elements);
                    }
                }
            }
            reader.close();
        /*    for (Map.Entry<String, Integer> entry : map.entrySet()) {
                    String key = entry.getKey();
                    Integer value = entry.getValue();
                    System.out.println("["+key+"]["+value+"]");
            }*/
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println("Bugger");
        }
        System.out.println("done");
    }

    public Instance newInstance(Instances dataType, String classType, String elements) {
        Instance instance;
        instance = new Instance(2);
        Attribute functionAtt = dataType.attribute("functions");
        instance.setValue(functionAtt, functionAtt.addStringValue(elements));
        instance.setDataset(dataType);
        instance.setClassValue(classType);
        return instance;
    }

    public void createTraingSet(String classType, int seed, double fraction) {
        ArrayList a = (ArrayList) hm.get(classType);
        Collections.shuffle(a, new Random(seed));
        int s = (int) (a.size() * fraction);
        for (int i = 0; i < s; i++) {
            Instance instance = this.newInstance(this.trainData, classType, (String) a.get(i));
            this.trainData.add(instance);
        }
        for (int i = s; i < a.size(); i++) {
            Instance instance = this.newInstance(this.testData, classType, (String) a.get(i));
            this.testData.add(instance);
        }
        ArrayList b = (ArrayList) excluded.get(classType);
        for (int i = 0; i < b.size(); i++) {
            Instance instance = this.newInstance(this.testData, classType, (String) b.get(i));
            this.testData.add(instance);
        }
        
    }

    public void loadDataset(String fileName) {
        try {
            BufferedReader reader = new BufferedReader(new FileReader(fileName));
            ArffReader arff = new ArffReader(reader);
            trainData = arff.getData();
            trainData.setClassIndex(1);
            reader.close();
        } catch (IOException e) {
            System.out.println("Problem found when reading: " + fileName);
        }
    }

    public void evaluate() {
        try {
            //for (int w = 10; w <= 10000; w *= 10) {
            StringToWordVector wordFilter = new StringToWordVector();
            wordFilter.setWordsToKeep(1000);
            //    System.out.println(w);
            WordTokenizer tokeniser = new WordTokenizer();
            String delimiters = " \t";
            tokeniser.setDelimiters(delimiters);
            wordFilter.setTokenizer(tokeniser);
            int[] atArray = {0};
            wordFilter.setAttributeIndicesArray(atArray);
            for (int c = 0; c < 2; c++) {
                classifier = new FilteredClassifier();
                classifier.setFilter(wordFilter);
                switch (c) {
                    case 0:
                        SMO smo = new SMO();
                        classifier.setClassifier(smo);
                        System.out.println("SMO");
                        break;
                    //   case 1:
                    //       NaiveBayes base = new NaiveBayes();
                    //       classifier.setClassifier(base);
                    //       System.out.println("NaiveBayes");
                    //       break;
                    //    case 2:
                    //        J48 j48 = new J48();
                    //       System.out.println("J48");
                    //       classifier.setClassifier(j48);
                    //      break;
                    case 1:
                        PART part = new PART();
                        classifier.setClassifier(part);
                        System.out.println("PART");
                        break;

                }
              //  for (int f = 10; f <= 90; f += 10) {
                    int reps = 10;
                    int f = 66;
                    //System.out.println((double) f / 100.0);
                    String[] classes = {"faculatative", "anaerobe", "aerobe"};
                    double[][] rocs = new double[classes.length][reps];
                    double[][][] matrix = new double[classes.length][classes.length][reps];
                    double[] pctCorr = new double[reps];
                    double[] wRocArea = new double[reps];
                    
                    for (int r = 0; r < reps; r++) {

                        trainData = initInstances(classes);
                        testData = initInstances(classes);
                        for (String s : classes) {
                            this.createTraingSet(s, r, (double) f / 100.0);
                        }
                        classifier.buildClassifier(trainData);
                        Evaluation eval = new Evaluation(testData);
                        eval.evaluateModel(classifier, testData);
                        for (int i = 0; i < eval.confusionMatrix().length; i++) {
                            rocs[i][r] = eval.areaUnderROC(i);
                            for (int j = 0; j < eval.confusionMatrix().length; j++) {
                                matrix[i][j][r] = eval.confusionMatrix()[i][j];
                            }
                        }
                        wRocArea[r] = eval.weightedAreaUnderROC();
                        pctCorr[r] = eval.pctCorrect();
                    }
                    System.out.println((double)trainData.numInstances()/(testData.numInstances()+trainData.numInstances())+"\t");
                        
                    for (int r = 0; r < reps; r++) {
                        for (int i = 0; i < classes.length; i++) {
                            for (int j = 0; j < classes.length; j++) {
                                System.out.print(matrix[i][j][r] + "\t");
                            }
                        }
                        for (int i = 0; i < classes.length; i++) {
                            System.out.print(rocs[i][r] + "\t");
                        }
                        System.out.print(wRocArea[r] + "\t");
                        System.out.println(pctCorr[r]);
                    }
                    System.out.println("\n\n\n");
               // }
            }
        } catch (Exception e) {
            System.out.println("Problem found when evaluating");
        }
    }

    public void classifierDetails() {
        System.out.println(classifier);
    }

    public void teachAndSave(String trainingFile) {
        System.out.println("Loading data");
        //this.loadDataset(trainingFile);
        this.loadRawData(trainingFile, "functions");
        this.evaluate();
        // this.classifierDetails();
    }

    public static void main(String[] args) {

        WekaClassifier learner = new WekaClassifier();
        //learner.teachAndSave("/home/ric/Programming/Classifiers/ClassifierPaper/Data/MetabolisamFunctionsFull.arff");
        learner.teachAndSave("/home/ric/Programming/Classifiers/ClassifierPaper/Data/classiferPaperFunctions010614.txt");

    }
}
