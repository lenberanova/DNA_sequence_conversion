#!/usr/bin/env python3

import hashlib
import os


def validate_input(input_file_path, input_md5, L):
    """check the path and file existence, number L value
    check the integrity of the file - comparison with expected MD5 value"""
    if not os.path.exists(input_file_path):
        print("The path '{}' does not exist.".format(input_file_path))
        return False

    if not os.path.isfile(input_file_path):
        print("'{}' is not a file.".format(input_file_path))
        return False

    with open(input_file_path, "rb") as input_file:
        content = input_file.read()
        encoded_content = hashlib.md5(content)

    if encoded_content.hexdigest() != input_md5:
        print("Invalid input file - integrity check failed.")
        return False

    if L > len(content) or L <= 0 or not isinstance(L, int):
        print("Number 'L' the value must be an integer in the range {}-{}.".format(1, len(content)))
        return False

    return True


def create_outup_element(fastq_element_list, byte_set):
    """create the sequence of basis and the sequence of quality score (length = L)"""
    dna_bases_dict = {0b00:"A", 0b01:"C", 0b10:"G", 0b11:"T"}
    sequence_of_bases = ""
    sequence_of_quality_score = ""

    for byte in byte_set:
        # bitwise right shift
        first_two_bits = byte >> 6
        base = dna_bases_dict[first_two_bits]
        sequence_of_bases = sequence_of_bases + base

        mask = 0b00111111
        # bitwise AND
        second_part = byte & mask
        quality_score = chr(second_part + 33)
        sequence_of_quality_score = sequence_of_quality_score + quality_score

    fastq_element = sequence_of_bases, sequence_of_quality_score
    fastq_element_list.append(fastq_element)
    return fastq_element_list


def process_input_file(input_file_path, L):
    """read the content of the binary file by number of bytes (L) and call the function that processes elements"""
    with open(input_file_path, "rb") as input_file:
        fastq_element_list = []
        byte_set = input_file.read(L)
        while len(byte_set) == L:
            create_outup_element(fastq_element_list, byte_set)
            byte_set = input_file.read(L)

        rest_count = len(byte_set)
        if rest_count > 0:
            create_outup_element(fastq_element_list, byte_set)
            print(rest_count)

        return fastq_element_list


def print_fastqt_data(output_fastq_data):
    """print the content of the input file on standard output in the FASTQ format"""
    for i in range(len(output_fastq_data)):
        print("@READ_{}".format(i + 1))
        print(output_fastq_data[i][0])
        print("+READ_{}".format(i + 1))
        print(output_fastq_data[i][1])


def execute_sequence_conversion(input_file_path, L):
    """execute sequence conversion and print on standard output"""
    fastq_data = process_input_file(input_file_path, L)
    print_fastqt_data(fastq_data)


input_path = "./dna_conversion_samples/input"
md5_input = "25a02f1331042be2856e652bda60e8de"
L = 7

is_valid = validate_input(input_path, md5_input, L)

if is_valid:
    execute_sequence_conversion(input_path, L)
