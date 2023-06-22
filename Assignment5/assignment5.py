#!/usr/local/bin/python3

import csv
import sys

import pyspark.sql.functions as f
from pyspark.sql import SparkSession
from pyspark.sql.window import Window


def question1(df):
    """
    Answering question 1:
    How many distinct protein annotations are found in the dataset?
    I.e. how many distinct InterPRO numbers are there?
    with pyspark
    :param df: pyspark dataframe
    :return: question number, answer, explanation
    """
    querry = df.select('_C11').distinct()
    result = querry.count()
    return 1, result, querry._jdf.queryExecution().toString()


def question2(df):
    """
    Answering question 2:
    How many annotations does a protein have on average?
    with pyspark
    :param df: pyspark dataframe
    :return: question number, answer, explanation
    """
    querry_grouped = df.groupby('_c0').count()
    querry = querry_grouped.agg({'count': 'mean'})
    result = querry.first()[0]
    return 2, result, querry._jdf.queryExecution().toString()


def question3(df):
    """
    Annswering question 3:
    What is the most common GO Term found?
    with pyspark
    :param df: pyspark dataframe
    :return: question number, answer, explanation
    """
    querry_grouped = df.where(df._c13 != '-').groupBy("_c13").count()
    querry_sorted = querry_grouped.sort("count", ascending=False)
    result = querry_sorted.first()[0]
    return 3, result, querry_sorted._jdf.queryExecution().toString()


def question4(df):
    """
    Answering question 4:
    What is the average size of an InterPRO feature found in the dataset?
    with pyspark
    :param df: pyspark dataframe
    :return: question number, answer, explanation
    """
    querry_grouped = df.groupby('_c11')
    querry_avg = querry_grouped.agg({'_c2' : 'mean'})
    result = querry_avg[0][0]
    return 4, result, querry_avg._jdf.queryExecution().toString()


def question5(df):
    """
    Answering question 5:
    What is the top 10 most common InterPRO features?
    with pyspark
    :param df: pyspark dataframe
    :return: question number, answer, explanation
    """
    querry_grouped = df.where(df._c11 != '-').groupby('_c11').count()
    querry_sorted = querry_grouped.sort('count', ascending=False)
    top_10 = [result[0] for result in querry_sorted.take(10)]
    return 5, top_10, querry_sorted._jdf.queryExecution().toString()


def question6(df):
    """
    Answering question 6:
    If you select InterPRO features that are almost the same size (within 90-100%)
    as the protein itself, what is the top10 then?
    with pyspark
    :param df: pyspark dataframe
    :return: question number, answer, explanation
    """
    querry = df.where(df._c11 != '-') \
        .withColumn('same_size', f.when((df['_c7'] - df['_c6']) > (df['_c2'] * 0.9), 1))
    querry_filtered = querry.filter(f.col('same_size').between(0, 2))
    querry_grouped = querry_filtered.groupBy('_c11').count()
    querry_sorted = querry_grouped.sort('count', ascending=False)
    top_10 = [result[0] for result in querry_sorted.take(10)]
    return 6, top_10, querry_sorted._jdf.queryExecution().toString()


def question7(df):
    """
    Answering question 7:
    If you look at those features which also have textual annotation,
    what is the top 10 most common word found in that annotation?
    with pyspark
    :param df: pyspark dataframe
    :return: question number, answer, explanation
    """
    querry_text = df.where(df._c12 != '-').withColumn('word', f.explode(f.split(f.col('_c12'), ' ')))
    querry_grouped = querry_text.groupBy('word').count()
    querry_sorted =  querry_grouped.sort('count', ascending=False)
    top_10 = [result[0] for result in querry_sorted.take(10)]
    return 7, top_10, querry_sorted._jdf.queryExecution().toString()


def question8(df):
    """
    Answering question 8:
    And the top 10 least common?
    with pyspark
    :param df: pyspark dataframe
    :return: question number, answer, explanation
    """
    querry_text = df.where(df._c12 != '-').withColumn('word', f.explode(f.split(f.col('_c12'), ' ')))
    querry_grouped = querry_text.groupBy('word').count()
    querry_sorted = querry_grouped.sort('count', ascending=True)
    top_10 = [result[0] for result in querry_sorted.take(10)]
    return 8, top_10, querry_sorted._jdf.queryExecution().toString()


def question9(df):
    """
    Answering question 9:
    Combining your answers for Q6 and Q7, what are the 10 most commons words found for the
    largest InterPRO features?
    with pyspark
    :param df: pyspark dataframe
    :return: question number, answer, explanation
    """
    querry = df.where(df._c11 != '-').withColumn('same_size', f.when((df['_c7'] - df['_c6']) > (df['_c2'] * 0.9), 1))
    querry_filtered = querry.filter(f.col('same_size').between(0, 2))
    querry_grouped = querry_filtered.withColumn('word', f.explode(f.split(f.col('_c12'), ' '))).groupBy('word')
    querry_sorted = querry_grouped.count().sort('count', ascending=False)
    top_10 = [result[0] for result in querry_sorted.take(10)]
    return 9, top_10, querry_sorted._jdf.queryExecution().toString()


def question10(df):
    """
    Answering question 10:
    What is the coefficient of correlation ($R^2$) between the size of the protein and
    the number of features found?
    with pyspark
    :param df: pyspark dataframe
    :return: question number, answer, explanation
    """
    querry_grouped = df.where(df._c11 != '-').select('_c0', '_c2', '_c11') \
        .withColumn('counts', f.count('_c11').over(Window.partitionBy('_c0')))
    querry_cleaned = querry_grouped.dropDuplicates(['_c0'])
    correlation_coefficient = querry_cleaned.stat.corr('_c2', 'counts')
    return 10, correlation_coefficient, querry_cleaned._jdf.queryExecution().simpleString()


def write_result(result):
    """
    Writes results of questions to an output csv file
    :param result: result object containing question number, answer
    and explanation
    """
    with open("output.csv", 'a', encoding='UTF-8', newline='') as output:
        csv_writer = csv.writer(output, dialect='excel')
        csv_writer.writerow(result)


def main(file):
    spark = SparkSession.builder.master('local[16]').appName('assignment5').getOrCreate()
    dataframe = spark.read.csv(file, sep='\t', header=False, inferSchema=True)
    write_result(question1(dataframe))
    write_result(question2(dataframe))
    write_result(question3(dataframe))
    write_result(question4(dataframe))
    write_result(question5(dataframe))
    write_result(question6(dataframe))
    write_result(question7(dataframe))
    write_result(question8(dataframe))
    write_result(question9(dataframe))
    write_result(question10(dataframe))
    spark.stop()


if __name__ == "__main__":
    main(sys.argv[1])
