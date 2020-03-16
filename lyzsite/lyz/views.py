import json
import traceback
import sys
import csv
import os
import re
from functools import reduce
from operator import and_
from django.shortcuts import render
from django import forms
from django.http import HttpResponse
from django.template import loader
from util_test import get_fasta


COLUMN_NAMES = dict(
    term='PDB ID',
    max_num='Maximum similar DNA',
    first ='picture',
    description = 'description',
)
NOPREF_STR = 'No preference'
RES_DIR = os.path.join(os.path.dirname(__file__), '..', 'res')


def _valid_result(res):
    '''
    Validate results returned by get_fasta.
    '''

    (HEADER, RESULTS) = [0, 1]
    ok = (isinstance(res, (tuple, list)) and
          len(res) == 2 and
          isinstance(res[HEADER], (tuple, list)) and
          isinstance(res[RESULTS], (tuple, list)))
    if not ok:
        return False

    n = len(res[HEADER])

    def _valid_row(row):
        return isinstance(row, (tuple, list)) and len(row) == n
    return reduce(and_, (_valid_row(x) for x in res[RESULTS]), True)


def _load_column(filename, col=0):
    """Load single column from csv file."""
    with open(filename) as f:
        col = list(zip(*csv.reader(f)))[0]
        return list(col)

def _load_res_column(filename, col=0):
    """Load column from resource directory."""
    return _load_column(os.path.join(RES_DIR, filename), col=col)


def _build_dropdown(options):
    """Convert a list to (value, caption) tuples."""
    return [(x, x) if x is not None else ('', NOPREF_STR) for x in options]

RANGE_WIDGET = forms.widgets.MultiWidget(widgets=(forms.widgets.NumberInput,
                                                  forms.widgets.NumberInput))


class SearchForm(forms.Form):

    def __init__(self, *args, **kwargs):
        super(SearchForm, self).__init__(
            *args, **kwargs)

    query = forms.CharField(
        label='Search Protein',
        help_text='e.g. 1J6Z',
        required=False)
    max_dna = forms.CharField(
    	label = 'Maximum DNA',
    	help_text='e.g. 20',
        required=False)
    label = forms.CharField('Maximum number of sequences',
        help_text='e.g. 20',
        required=False)
    e_value = forms.CharField(
        label = 'E value',
        help_text='e.g. evalue',
        required=False)    
    show_args = forms.BooleanField(label='Show args_to_ui',
                                   required=False)



def index(request):
    context = {}
    template = loader.get_template('index.html')
    res = None
    if request.method == 'GET':     
        form = SearchForm(request.GET)
        if form.is_valid():
            args = {}
            print(form.cleaned_data['query'])
            if form.cleaned_data['query']:
                args['terms'] = form.cleaned_data['query']
            print(args)
            if form.cleaned_data['max_dna']:
                args['Maximum DNA'] = form.cleaned_data['max_dna']
                args['Maximum number of sequences'] = form.cleaned_data['max_dna']
            if form.cleaned_data['e_value']:
                args['E value'] = form.cleaned_data['e_value']
            if form.cleaned_data['show_args']:
                context['args'] = 'args_to_ui = ' + json.dumps(args, indent=2)
            try:
                res = get_fasta(args)
            except Exception as e:
                print('Exception caught')
                bt = traceback.format_exception(*sys.exc_info()[:3])
                context['err'] = """
                An exception was thrown in:
                <pre>{}
{}</pre>
                """.format(e, '\n'.join(bt))

                res = None
    else:
        form = SearchForm()
    if res is None:
        context['result'] = None
    else:
        columns, result = res
        if result and isinstance(result[0], str):
            result = [(r,) for r in result]
        context['result'] = result
        context['num_results'] = len(result)
        context['form'] = form
    return render(request, 'index.html', context)
