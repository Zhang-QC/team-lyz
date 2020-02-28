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


# Create your views here.

#def index(request):
    #return HttpResponse("Hello, world. You're at the lyz index.")
    #return SearchForm(request)
def _valid_result(pdb):
    '''
    Validate results returned by find_courses.

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
    '''

    if re.search("[A-Z0-9]{4}",res) != None and len(res) == 4:
        return True
    return False



class SearchForm(forms.Form):
    def __init__(self, *args, **kwargs):
        super(SearchForm, self).__init__(
            *args, **kwargs)


    query = forms.CharField(
        label='Search Protein',
        help_text='e.g. 1J6Z',
        required=True)

def index(request):
    context = {}
    template = loader.get_template('index.html')
    res = None
    if request.method == 'GET':
        # create a form instance and populate it with data from the request:
        
        print(request.GET)
        if "terms" not in request.GET:
            print("skip")
            #d = template.render(request, context)
            #print(d)
            print('return suuc')
            return HttpResponse(template.render({}, request)) #template.render(context, request)
        print("NotSKIP")
        
        form = SearchForm(request.GET)
        print("FORM_test")
        # check whether it's valid:
        if form.is_valid():
            # Convert form data to an args dictionary for find_courses
            if 'query' not in form.cleaned_data:
                print("skip")
                return HttpResponse(template.render({}, request))
            print("Not skip")
            args = {}
            if form.cleaned_data['query']:
                args['terms'] = form.cleaned_data['query']
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

    # Handle different responses of res
    if res is None:
        context['result'] = None
    elif isinstance(res, str):
        context['result'] = None
        context['err'] = res
        result = None
    elif not _valid_result(res):
        context['result'] = None
        context['err'] = ('Return of XXX has the wrong data type.')
    else:
        columns, result = res

        # Wrap in tuple if result is not already
        if result and isinstance(result[0], str):
            result = [(r,) for r in result]

        context['result'] = result
        context['num_results'] = len(result)
        #context['columns'] = [COLUMN_NAMES.get(col, col) for col in columns]

    context['form'] = form
    #return render(request, 'index.html', context)
    return HttpResponse(template.render({}, request))
