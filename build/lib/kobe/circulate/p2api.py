"""
This file provides an example binding to ESO's RESTful phase 2 programmatic interface.
The data required and returned by many of the API calls is instrument-specific.
In order to fully understand the data and behaviour of the API please consult
the documentation at https://www.eso.org/copdemo/apidoc/
The api URL specified with each API should help to locate the corresponding
API documentation at the above URL.

Unless otherwise mentioned, each API call consistently returns a tuple (data, version):
    * data    - returned data in JSON format
    * version - version of the data, required for future modifications (HTTP ETags are used for this versioning)

___________________________________________________________________________
MIT License
___________________________________________________________________________
Copyright (c) 2018 Thomas Bierwirth, European Southern Observatory

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
from __future__ import print_function
import os
import json
import requests

# $Id: p2api.py 235106 2019-01-23 15:09:20Z tbierwir $
P2_API_VERSION = '0.11'

API_URL = {
    'production': 'https://www.eso.org/cop/api/v1',
    'demo': 'https://www.eso.org/copdemo/api/v1',
}
LOGIN_URL = {
    'production': 'https://www.eso.org/cop/api/login',
    'demo': 'https://www.eso.org/copdemo/api/login',    
}
P2FC_URL = {
    'production': 'https://www.eso.org/p2fc/',
    'demo': 'https://www.eso.org/p2fcdemo/',
}

class P2Error(Exception):
    pass

class ApiConnection(object):
# ---------- PUBLIC API ----------
    def __init__(self, environment, username, password, debug=False):
        """
        Set environment, log in to the API and obtain access token for further
        calls. Returns authenticated api connection required to invoke other
        API calls.

        usage:
           ::

              api = p2api.ApiConnection('production', 52052', 'tutorial')
        """
        self.debug = debug
        self.request_count = 0
        self.apiUrl = API_URL[environment]
        self.loginUrl = LOGIN_URL[environment]    
        self.p2fc_url = P2FC_URL[environment]
        r = requests.post(
            self.loginUrl,
            data={'username': username, 'password': password})
        if r.status_code == requests.codes.ok:
            body = r.json()
            self.access_token = body['access_token']
        else:
            raise P2Error(r.status_code, 'POST', self.loginUrl, 'cannot login')

# ---------- RUN-LEVEL APIs ----------
    def getRuns(self):
        """Retrieve observing runs owned by or delegated to logged-in user.

        api: GET /obsRuns

        usage:
           ::

              runs, _ = api.getRuns()
        """
        return self.get('/obsRuns')

    def getRun(self, runId):
        """Get a single observing run.

        api: GET /obsRuns/{runId}

        usage:
           ::

              run, _ = api.getRun(runId)
        """
        return self.get('/obsRuns/%d' % runId)

    def getPhase1Targets(self, runId):
        """Retrieve the list of targets defined in phase 1.

        api: GET /obsRuns/{runId}/phase1/targets

        usage:
           ::

              targets, _ = api.getPhase1Targets(runId)
        """
        return self.get('/obsRuns/%d/phase1/targets' % runId)

    def getBlockedTimes(self, runId):
        """Retrieve a list of blocked times on the telescope in the
        scheduled period of this observing run.

        api: GET /obsRuns/{runId}/phase1/blockedTimes

        usage:
           ::

              targets, _ = api.getBlockedTimes(runId)
        """
        return self.get('/obsRuns/%d/phase1/blockedTimes' % runId)

    def submitRun(self, runId):
        """Submit the observing run to ESO for review.

        api  : POST /obsRuns/{runId}/submit

        usage:
           ::

              api.submitRun(runId)
        """
        return self.post('/obsRuns/%d/submit' % runId)

    def getReadmeSchema(self, runId):
        """Retrieve the schema definition for the ReadMe file of this observing run.

        api: GET /obsRuns/{runId}/readme/schema

        usage:
           ::

              readmeSchema, _ = api.getReadmeSchema(runId)
        """
        return self.get('/obsRuns/%d/readme/schema' % runId)

    def getReadme(self, runId):
        """Retrieve the ReadMe file of this observing run.

        api: GET /obsRuns/{runId}/readme

        usage:
           ::

              readme, readmeVersion = api.getReadme(runId)
        """
        return self.get('/obsRuns/%d/readme' % runId)

    def saveReadme(self, runId, readme, version):
        """
        Update the ReadMe file.

        api: PUT /obsRuns/{runId}/readme.

        usage:
           ::

              readme, readmeVersion = api.saveReadme(runId, readme, version)
        """
        return self.put('/obsRuns/%d/readme' % runId, readme, version)

# ---------- CONTAINER-LEVEL APIs ----------
    def getContainer(self, containerId):
        """Get container information.

        api: GET /containers/{containerId}

        usage:
           ::

              container, containerVersion = api.getContainer(containerId)
        """
        return self.get('/containers/%d' % containerId)

    def saveContainer(self, container, version):
        """Update container information.

        api: PUT /containers/{containerId}

        usage:
           ::

              container, containerVersion = api.getContainer(containerId)
              # modify some container properties ...
              container['userPriority'] = 10
              # save changes
              container, containerVersion = api.saveContainer(container, containerVersion)
        """
        return self.put('/containers/%d' % container['containerId'], container, version)

    def deleteContainer(self, containerId, version):
        """Delete container. The container must be empty.

        api: DELETE /containers/{containerId}

        usage:
           ::

              container, containerVersion = api.getContainer(container)
              api.deleteContainer(container['containerId'], containerVersion)
        """
        return self.delete('/containers/%d' % containerId, etag=version)

    def createItem(self, itemType, containerId, name):
        """Create a new OB, CB, Folder, Group, Concatenation or TimeLink
        at the end of the given container.

        The initial values for the observing constraints of a new OB are those
        specified during phase 1, or the instrument defaults.

        Visitor Mode: Groups, Concatenations and TimeLinks cannot be created.

        Service Mode: Groups, Concatenations and TimeLinks cannot contain other
        containers, only OBs and CBs.

        See also wrappers below to create specific items by name.

        api: POST /containers/{containerId}/items

        usage:
           ::

              ob, obVersion = api.createItem('OB', containerId, 'OliOB')
        """
        return self.post('/containers/%d/items' % containerId, {'itemType': itemType, 'name': name})

    def createOB(self, containerId, name):
        """Create a new Observing Block.

        api: POST /containers/{containerId}/items

        usage:
           ::

              ob, obVersion = api.createOB(containerId, 'OliOB')
        """
        return self.createItem('OB', containerId, name)

    def createCB(self, containerId, name):
        """Create a new Calibration Block.

        api  : POST /containers/{containerId}/items

        usage:
           ::

              cb, cbVersion = api.createCB(containerId, 'OliCB')
        """
        return self.createItem('CB', containerId, name)

    def createFolder(self, containerId, name):
        """Create a new Folder.

        api: POST /containers/{containerId}/items

        usage:
           ::

              fld, fldVersion = api.createFolder(containerId, 'OliFolder')
        """
        return self.createItem('Folder', containerId, name)

    def createGroup(self, containerId, name):
        """Create a new Group.

        api: POST /containers/{containerId}/items

        usage:
           ::

              grp, grpVersion = api.createGroup(containerId, 'OliGroup')
        """
        return self.createItem('Group', containerId, name)

    def createConcatenation(self, containerId, name):
        """Create a new Concatenation.

        api: POST /containers/{containerId}/items

        usage:
           ::

              con, conVersion = api.createConcatenation(containerId, 'OliConcatenation')
        """
        return self.createItem('Concatenation', containerId, name)

    def createTimeLink(self, containerId, name):
        """Create a new TimeLink

        api: POST /containers/{containerId}/items

        usage:
           ::

              tl, tlVersion = api.createTimeLink(containerId, 'OliTimeLink')
        """
        return self.createItem('TimeLink', containerId, name)

    def getItems(self, containerId):
        """ Get ordered list of OBs and containers inside a container.

        api: GET /containers/{containerId}/items

        usage:
           ::

              items, itemsVersion = api.getItems(containerId)
        """
        return self.get('/containers/%d/items' % containerId)

    def reorderItems(self, containerId, items, version):
        """Change order of OBs and containers inside a container.

        api: PUT /containers/{containerId}/items

        usage:
           ::

              items, itemsVersion = api.getItems(containerId)
              # modify items order ...
              # save changes
              items, itemsVersion = api.reorderItems(containerId, items, itemsVersion)
        """
        return self.put('/containers/%d/items' % containerId, items, version)

    def moveFolder(self, folderId, containerId):
        """Move existing folder into destination container.

        api: POST /containers/{containerId}/items/append

        usage:
           ::

              api.moveFolder(folderId, containerId)
        """
        return self.post('/containers/%d/items/append' % containerId,
                         {'itemType': 'Folder', 'containerId': folderId})

    def moveGroup(self, groupId, containerId):
        """Move existing group into destination container.

        api: POST /containers/{containerId}/items/append

        usage:
           ::

              api.moveGroup(groupId, containerId)
        """
        return self.post('/containers/%d/items/append' % containerId,
                         {'itemType': 'Group', 'containerId': groupId})

    def moveConcatenation(self, concatenationId, containerId):
        """Move existing concatenation into destination container.

        api: POST /containers/{containerId}/items/append

        usage:
           ::

              api.moveConcatenation(concatenationId, containerId)
        """
        return self.post('/containers/%d/items/append' % containerId,
                         {'itemType': 'Concatenation', 'containerId': concatenationId})

    def moveTimeLink(self, timeLinkId, containerId):
        """Move existing timelink into destination container.

        api: POST /containers/{containerId}/items/append

        usage:
           ::

              api.moveTimeLink(folderId, containerId)
        """
        return self.post('/containers/%d/items/append' % containerId,
                         {'itemType': 'TimeLink', 'containerId': timeLinkId})

    def moveOB(self, obId, containerId):
        """Move existing OB into destination container.

        api: POST /containers/{containerId}/items/append

        usage:
           ::

              api.moveOB(obId, containerId)
        """
        return self.post('/containers/%d/items/append' % containerId,
                         {'itemType': 'OB', 'obId': obId})

    def moveCB(self, obId, containerId):
        """Move existing CB into destination container.

        api: POST /containers/{containerId}/items/append

        usage:
           ::

              api.moveCB(obId, containerId)
        """
        return self.post('/containers/%d/items/append' % containerId,
                         {'itemType': 'CB', 'obId': obId})

    def getOBSummary(self, containerId):
        """Get more detailed summary of all OBs in a container

        api: GET /containers/{containerId}/items/summary

        usage:
           ::

              api.getOBSummary(containerId)
        """
        return self.get('/containers/%d/items/summary' % containerId)

    def verifyContainer(self, containerId, submit):
        """Verify a Group, Concatenation or TimeLink.

        If submit is false, a Group, Concatenation, or TimeLink can be verified
        in any container status. If submit is true, the container must have
        container status "-" or "P".

        OBs inside the container must be verified first, i.e. they have a status
        of "D" or later.

        The instrument-specific verification scripts for a container are executed.
        If submit is true and verification is successful, the container status
        is changed to "D".

        A change of the container status will invalidate the ETag of the container.

        api: POST /containers/{containerId}/verify

        usage:
           ::

              response, _ = api.verifyContainer(containerId, True)
              if response['observable']:
                 print('Your container is observable!')
              else:
                 print('Your container is >>not observable<<. See messages below.')
                 print(' ')
                 print(' '.join(response['messages']))
        """
        return self.post('/containers/%d/verify' % containerId, {'submit': submit})

    def reviseContainer(self, containerId):
        """Open Group/Concatenation/TimeLink for further changes.

        Change the container status back to "P".
        This is the counterpart to verification with submission.

        api: POST /containers/{containerId}/revise

        usage:
           ::

              api.reviseContainer(containerId)
        """
        return self.post('/containers/%d/revise' % containerId)

    def exportContainer(self, containerId, filename):
        """Export this container into a ZIP file.

        api: GET /containers/{containerId}/export

        usage:
           ::

              api.exportContainer(containerId, 'myContainer.zip')
        """
        self.downloadBinaryFile('/containers/%d/export' % containerId, filename)

# ---------- OB-LEVEL APIs ----------
    def getOB(self, obId):
        """Get existing OB or CB.

        api: GET /obsBlocks/{obId}

        usage:
           ::

              ob, obVersion = api.getOB(obId)
        """
        return self.get('/obsBlocks/%d' % obId)

    def saveOB(self, ob, version):
        """Update existing OB or CB.

        The format of an OB is instrument-specific and differs for OBs and CBs.
        The instrument defines the applicable observingConstraints.
        A CB does not have a target. The properties itemType, obId, obStatus,
        ipVersion, exposureTime and executionTime cannot be changed.
        The OB must have the latest ipVersion. Refer to the API spec
        for further updating constraints.

        api: PUT /obsBlocks/{obId}

        usage:
           ::

              ob, obVersion = api.getOB(obId)
              # modify some OB properties ...
              ob['userPriority'] = 42
              ob['target']['name']         = 'M 32 -- Interacting Galaxies'
              ob['target']['ra']           = '00:42:41.825'
              ob['target']['dec']          = '-60:51:54.610'
              ob['constraints']['name']    = 'My hardest constraints ever'
              ob['constraints']['airmass'] = 1.3
              # save changes
              ob, obVersion = api.saveOB(ob, obVersion)
        """
        return self.put('/obsBlocks/%d' % ob['obId'], ob, version)

    def deleteOB(self, obId, version):
        """Delete existing OB or CB.

        Visitor Mode: The OB must have never been executed.

        Service Mode: The OB must have status (-), (P) or (D) and must have never been executed.

        api: DELETE /obsBlocks/{obId}

        usage:
           ::

              ob, obVersion = api.getOB(obId)
              api.deleteOB(ob, obVersion)
        """
        return self.delete('/obsBlocks/%d' % obId, etag=version)

    def verifyOB(self, obId, submit):
        """Verify OB and submit it for observation.

        If submit is false, an OB can be verified in any status. If submit is
        true, a loose service mode OB must be in status "P" or "-" and a visitor
        mode OB must be in a status other than "S", "A", "X", or "+".

        api: POST /obsBlocks/{obId}/verify

        usage:
           ::

              response, _ = api.verifyOB(obId, True)
              if response['observable']:
                 print('Your OB' , obId, ob['name'], 'is observable!')
              else:
                 print('OB', obId, 'is >>not observable<<. See messages below.')
                 print(' ')
                 print(' '.join(response['messages']))
        """
        return self.post('/obsBlocks/%d/verify' % obId, {'submit': submit})

    def reviseOB(self, obId):
        """Back to the drawing board.

        Change the OB status back to "P", i.e. mark it as not observable.
        This is the counterpart to verification with submission.

        A visitor mode OB can be in any status. The status is changed to "P".
        The OB will be removed from any visitor execution sequences because
        it is no longer in an observable status; this changes the ETag of
        the affected visitor execution sequences.

        A service mode OB can only be revised in status "D". The status is
        changed to "P". If the OB is a member of a Group, Concatenation or
        TimeLink with container status other than "-" or "P", the container
        status will be changed to "P".

        This call invalidates the ETag of the OB, and possibly of the parent
        container, if the container status was changed.

        api: POST /obsBlocks/{obId}/revise

        usage:
           ::

              pi.reviseOB(obId)
        """
        return self.post('/obsBlocks/%d/revise' % obId)

    def duplicateOB(self, obId, containerId=0):
        """Duplicate OB.

        Create a new OB or CB identical to this one in the same or a different
        container. The new OB is appended at the end of the container. The
        destination may be in a different observing run of the same instrument.

        api: POST /obsBlocks/{obId}/duplicate

        usage:
           ::

              duplicatedOB, duplicatedOBVersion = api.duplicateOB(obId)
              duplicatedOB, duplicatedOBVersion = api.duplicateOB(obId, containerId)
        """
        data = {}
        if containerId > 0:
            data = {'containerId': containerId}
        return self.post('/obsBlocks/%d/duplicate' % obId, data)

    def migrateOB(self, obId, dryRun):
        """Migrate OB to latest instrument package.

        api: POST /obsBlocks/{obId}/migrate

        usage:
           ::

              result, _ = api.migrateOB(obId, false)
              result, _ = api.migrateOB(obId, true)
        """
        return self.post('/obsBlocks/%d/migrate' % obId, {'dryRun': dryRun})

    def exportOB(self, obId, filename):
        """Export this OB into a ZIP file.

        api: GET /obsBlocks/{obId}/export

        usage:
           ::

              api.exportOB(obId, 'myOB.zip')
        """
        self.downloadBinaryFile('/obsBlocks/%d/export' % obId, filename)

    def getOBExecutions(self, obId, night):
        """Get executions during a given night.

        api: GET /obsBlocks/{obId}/executions

        usage:
           ::

              obExecutions, _ = getOBExecutions(obId, '2018-01-01')
        """
        return self.get('/obsBlocks/%d/executions?night=%s' % (obId, night))

    def getOBCompletelyFilledIn(self, obId):
        """Check that the OB's templates are completely filled in.

        Performs a quick check that there are neither missing templates
        nor missing parameters on this OB.

        This is not a full validation.

        api: GET /obsBlocks/{obId}/templates/completelyFilledIn

        usage:
           ::

              response, _ = api.getOBCompletelyFilledIn(obId)
              print('OB complete:', response['complete'])
        """
        return self.get('/obsBlocks/%d/templates/completelyFilledIn' % obId)

    def getAbsoluteTimeConstraints(self, obId):
        """Get list of absolute time constraints.

        The list is sorted in ascending order by from.
        Date and time are in ISO-8601 format for the UTC timezone.
        A CB does not have absolute time constraints and returns an empty list.

        api: GET /obsBlocks/{obId}/timeConstraints/absolute

        usage:
           ::

              absTCs, atcVersion = api.getAbsoluteTimeConstraints(obId)
        """
        return self.get('/obsBlocks/%d/timeConstraints/absolute' % obId)

    def saveAbsoluteTimeConstraints(self, obId, timeConstraints, version):
        """Set list of absolute time constraints.

        Date and time are in ISO-8601 format and must use UTC timezone.
        Time intervals must not overlap and from must be earlier than to.

        api: PUT /obsBlocks/{obId}/timeConstraints/absolute

        usage:
           ::

              absTCs, atcVersion = api.getAbsoluteTimeConstraints(obId)
              api.saveAbsoluteTimeConstraints(obId,[
                {
                    'from': '2017-09-01T00:00',
                    'to': '2017-09-30T23:59'
                },
                {
                    'from': '2017-11-01T00:00',
                    'to': '2017-11-30T23:59'
                }
              ], atcVersion)
        """
        return self.put('/obsBlocks/%d/timeConstraints/absolute' % obId, timeConstraints, version)

    def getSiderealTimeConstraints(self, obId):
        """Get list of sidereal time constraints.

        The list is sorted in ascending order by from.
        A CB does not have sidereal time constraints and returns an empty list.

        api: GET /obsBlocks/{obId}/timeConstraints/sidereal

        usage:
           ::

              sidTCs, stcVersion = api.getSiderealTimeConstraints(obId)
        """
        return self.get('/obsBlocks/%d/timeConstraints/sidereal' % obId)

    def saveSiderealTimeConstraints(self, obId, timeConstraints, version):
        """Set list of sidereal time constraints.

        Time intervals must not overlap. A time interval ending at the starting time
        of another is considered an overlap; avoid this by combining the intervals
        into a single one.

        There can be a single time interval ending at 24:00; it will not overlap
        a time interval starting at 0:00. An interval ending at 24:00 and one
        starting at 00:00 will be combined into a single interval.

        24:00 cannot be used as starting time, and 00:00 cannot be used
        as ending time.

        api: PUT /obsBlocks/{obId}/timeConstraints/sidereal

        usage:
           ::

              sidTCs, stcVersion = api.getSiderealTimeConstraints(obId)
              api.saveSiderealTimeConstraints(obId,[
                {
                    'from': '00:00',
                    'to': '01:00'
                },
                {
                    'from': '03:00',
                    'to': '05:00'
                }
              ], stcVersion)
        """
        return self.put('/obsBlocks/%d/timeConstraints/sidereal' % obId, timeConstraints, version)

    def createTemplate(self, obId, name):
        """Attach a new template.

        An OB can have only one acquisition template. Attempting to add a
        second one will result in an error.

        The OB must have the latest ipVersion.

        Visitor Mode: The OB can be changed in any status.

        Service Mode: The OB must have status (-) or (P).

        api: POST /obsBlocks/{obId}/templates

        usage:
           ::

              tpl, tplVersion = api.createTemplate(obId, 'UVES_blue_acq_slit')
        """
        return self.post('/obsBlocks/%d/templates' % obId, {'templateName': name})

    def getTemplates(self, obId):
        """Get all attached templates.

        api: GET /obsBlocks/{obId}/templates

        usage:
           ::

              templates, templatesVersion = api.getTemplates(obId)
              for t in templates:
                 print ('template found:', t['templateName'])
        """
        return self.get('/obsBlocks/%d/templates' % obId)

    def reorderTemplates(self, obId, templates, version):
        """Change order of attached templates.

        The list must contain the same items as in the previous GET request.
        Only the order can be changed. Templates must be explicitly added
        with a POST request or deleted with a DELETE request. If there is
        an acquisition template, it must be the first in the list.

        The OB must have the latest ipVersion.

        Visitor Mode: The OB can be changed in any status.

        Service Mode: The OB must have status (-) or (P).

        api: PUT /obsBlocks/{obId}/templates

        usage:
           ::

              tpls, tplsVersion = api.getTemplates(obId)
              # change order of tpls ...
              tpls, tplsVersion = api.reorderTemplates(obId, tpls, tplsVersion)
        """
        return self.put('/obsBlocks/%d/templates' % obId, templates, version)

    def getTemplate(self, obId, templateId):
        """Get existing template.

        A parameter value may be null if there is no default value and a
        value was never set.

        For parameters of type 'paramfile' or 'file' a string is returned.
        For a 'paramfile' this is the value of PAF.NAME in the PAF header;
        for a 'file' the first few characters are returned.
        File parameters are set with a separate call.

        api: GET /obsBlocks/{obId}/templates/{templateId}

        usage:
           ::

              tpl, tplVersion = api.getTemplate(obId, templateId)
        """
        return self.get('/obsBlocks/%d/templates/%d' % (obId, templateId))

    def setTemplateParams(self, obId, template, params, version):
        """Helper API to make it easier to update template parameters.
           Uses saveTemplate() under the hood.


        usage:
           ::

              tpl, tplVersion = api.createTemplate(obId, 'UVES_blue_acq_slit')
              tpl, tplVersion  = api.setTemplateParams(obId, tpl, {
                'TEL.GS1.ALPHA': '11:22:33.000',
                'INS.COLL.NAID': 'COLL_SR+7',
                'INS.DPOL.MODE': 'ON'
              }, tplVersion)
        """
        for p in template['parameters']:
            p['value'] = params.get(p['name'], p['value'])
        return self.saveTemplate(obId, template, version)

    def saveTemplate(self, obId, template, version):
        """Update existing template.

        The templateId, templateName and type cannot be changed.

        When GET returns a value of null, it can be changed to a proper value
        (or left at null). But a proper value cannot be unset to null.

        File parameters must be set with a separate call; the string value
        returned from GET cannot be changed.

        The OB must have the latest ipVersion.

        Visitor Mode: The OB can be changed in any status.

        Service Mode: The OB must have status (-) or (P).

        api: PUT /obsBlocks/{obId}/templates/{templateId}

        usage:
           ::

              It is recommended to use setTemplateParams() instead
        """
        return self.put('/obsBlocks/%d/templates/%d' % (obId, template['templateId']), template, etag=version)

    def getFileParam(self, obId, templateId, param, filename):
        """Get value of a template parameter of type "file" or "paramfile".

        api: GET /obsBlocks/{obId}/templates/{templateId}/{fileParameter}

        usage:
           ::

              _, version = api.getFileParam(obId, templateId, 'SEQ.REF.FILE1', myFile.paf')
        """
        return self.downloadTextFile('/obsBlocks/%d/templates/%d/%s' % (obId, templateId, param), filename)

    def saveFileParam(self, obId, templateId, param, filename, version):
        """Set value of a parameter of type "file" or "paramfile".

        api: PUT /obsBlocks/{obId}/templates/{templateId}/{fileParameter}

        usage:
           ::

              _, version = api.getFileParam(obId, templateId, 'SEQ.REF.FILE1', myFile.paf')
              api.saveFileParam(obId, templateId, 'SEQ.REF.FILE1', newFile.paf', version)
        """
        return self.uploadFile('PUT', '/obsBlocks/%d/templates/%d/%s' % (obId, templateId, param), filename,
                               'text/plain', version)

    def deleteTemplate(self, obId, templateId, version):
        """Delete existing template.

        The OB must have the latest ipVersion.

        Visitor Mode: The OB can be changed in any status.

        Service Mode: The OB must have status (-) or (P).

        api: DELETE /obsBlocks/{obId}/templates/{templateId}

        usage:
           ::

              _, tplVersion = api.getTemplate(obId, templateId)
              api.deleteTemplate(obId, templateId, tplVersion)
        """
        return self.delete('/obsBlocks/%d/templates/%d' % (obId, templateId), etag=version)

    def duplicateTemplate(self, obId, templateId, destinationObId=0):
        """Duplicate this template.

        Create a new template identical to this one in this or a different OB.
        The new template is appended at the end of the target OB's template list.

        The new template will have a new templateId. All other fields, including
        any file and paramfile parameters, will be the same as in the original template.

        An OB can have only one acquisition template. Attempting to add a second
        one will result in an error.

        api: POST /obsBlocks/{obId}/templates/{templateId}/duplicate

        usage:
           ::

              duplicatedTemplate, duplicatedTemplateVersion = api.duplicateTemplate(obId, templateId)
              duplicatedTemplate, duplicatedTemplateVersion = api.duplicateTemplate(obId, templateId, destinationObId)
        """
        data = {}
        if destinationObId > 0:
            data = {'obId': destinationObId}
        return self.post('/obsBlocks/%d/templates/%d/duplicate' % (obId, templateId), data)

    def getEphemerisFile(self, obId, filename):
        """Get the ephemeris file of this OB and save it into the given filename.

        If the ephemeris file was never saved, or was deleted, an empty file
        is returned.

        api: GET /obsBlocks/{obId}/ephemeris

        usage:
           ::

              _, ephVersion = api.getEphemerisFile(obId, 'ephem.txt')
        """
        url = '/obsBlocks/%d/ephemeris' % obId
        return self.downloadTextFile(url, filename)

    def saveEphemerisFile(self, obId, filename, version):
        """Save the ephemeris file of this OB.

        A CB cannot have an ephemeris file.

        Service Mode: The OB must have status (-) or (P).

        api: PUT /obsBlocks/{obId}/ephemeris

        usage:
           ::

              _, ephVersion = api.getEphemerisFile(obId, 'delete_me.txt')
              _, ephVersion = api.saveEphemerisFile(obId, 'ephem.txt', ephVersion)
        """
        return self.uploadFile('PUT', '/obsBlocks/%d/ephemeris' % obId, filename,
                               'text/plain', version)

    def deleteEphemerisFile(self, obId, version):
        """Delete the ephemeris file of this OB.

        Service Mode: The OB must have status (-) or (P).

        api: DELETE /obsBlocks/{obId}/ephemeris

        usage:
           ::

              _, ephVersion = api.getEphemerisFile(obId, 'delete_me.txt')
              api.deleteEphemerisFile(obId, ephVersion)
        """
        return self.request('DELETE', '/obsBlocks/%d/ephemeris' % obId, etag=version)

    def getFindingChartNames(self, obId):
        """List attached finding charts.

        Returns the list of names of the attached finding charts in the order
        of their index. If no finding charts are attached, an empty list will
        be returned.

        api: GET /obsBlocks/{obId}/findingCharts

        usage:
           ::

              fcNames, _ = api.getFindingChartNames(obId)
        """
        return self.get('/obsBlocks/%d/findingCharts' % obId)

    def addFindingChart(self, obId, filename):
        """Attach a new finding chart in format image/jpeg.

        Up to 5 finding charts can be attached to a single OB.
        The maximum file size is 1MB (2^20 bytes).

        The OB must have the latest ipVersion.

        Service Mode: The OB must have status (-) or (P).

        api: POST /obsBlocks/{obId}/findingCharts

        usage:
           ::

              api.addFindingChart(obId, filename)
        """
        return self.uploadFile('POST', '/obsBlocks/%d/findingCharts' % obId, filename,
                               'image/jpeg')

    def generateFindingChart(self, obId):
        """Generate and new finding chart(s) in format image/jpeg using ESO's p2fc service.

        Up to 5 finding charts can be attached to a single OB.
        The maximum file size is 1MB (2^20 bytes).

        The OB must have the latest ipVersion.

        Service Mode: The OB must have status (-) or (P).

        api: POST /obsBlocks/{obId}/findingCharts

        usage:
           ::

              api.generateFindingChart(obId)
        """
        data = {
            'url': self.apiUrl + '/obsBlocks/%d' % obId,
            'access_token': self.access_token
        }
        r = requests.post(
            self.p2fc_url,
            data = data
        )
        content_type = r.headers['Content-Type'].split(';')[0]
        etag = r.headers.get('ETag', None)
        # handle response
        if 200 <= r.status_code < 300:
            if content_type == 'application/json':
                data = r.json()
                return data, etag
            return None, etag
        elif content_type == 'application/json' and 'error' in r.json():
            raise P2Error(r.status_code, 'POST', self.p2fc_url, r.json()['error'])
        else:
            raise P2Error(r.status_code, 'POST', self.p2fc_url, r.text)
        '''
        return self.uploadFile('POST', '/obsBlocks/%d/findingCharts' % obId, filename,
                               'image/jpeg')
        '''
    def getFindingChart(self, obId, index, filename):
        """Download finding chart save it into the given filename.

        api: GET /obsBlocks/{obId}/findingCharts/{index}

        usage:
           ::

              api.getFindingChart(obId, 1, 'fc1.jpg')
        """
        url = '/obsBlocks/%d/findingCharts/%d' % (obId, index)
        return self.downloadBinaryFile(url, filename)

    def deleteFindingChart(self, obId, index):
        """Delete finding chart.

        The finding chart at the given index is deleted.
        Finding charts with a higher index 'slide down',
        i.e. finding charts are always numbered consecutively
        starting with 1.

        The OB must have the latest ipVersion.

        Service Mode: The OB must have status (-) or (P).

        api: DELETE /obsBlocks/{obId}/findingCharts/{index}

        usage:
           ::

              api.deleteFindingChart(obId, 1)
        """
        return self.delete('/obsBlocks/%d/findingCharts/%d' % (obId, index))

    def duplicateFindingChart(self, obId, index, destinationObId):
        """Duplicate this finding chart to another OB.

       Attach this finding chart to a different OB.
       The finding chart is appended at the end of the target OB's
       finding chart list.

        The target OB must have the latest ipVersion.

        Visitor Mode: The target OB can have any status.

        Service Mode: The target OB must have status (-) or (P).

        api: POST /obsBlocks/{obId}/findingCharts/{index}/duplicate

        usage:
           ::

              api.duplicateFindingChart(obId, 1, destinationObId)
        """
        return self.post('/obsBlocks/%d/findingCharts/%d/duplicate' % (obId, index), {'obId': destinationObId})

    def getOBSchedulingInfo(self, obId):
        """Get OB scheduling information.

        api: GET /obsBlocks/{obId}/schedule

        usage:
           ::

               obInfo, obInfoVersion = api.getOBSchedulingInfo(obId)
        """
        return self.get('/obsBlocks/%d/schedule' % obId)

    def saveOBSchedulingInfo(self, obId, obInfo, version):
        """Update OB scheduling information.

        api: PUT /obsBlocks/{obId}/schedule

        usage:
           ::

               obInfo, obInfoVersion = api.getOBSchedulingInfo(obId)
               # change something
               obInfo['targetName'] = 'My Best Target'
               api.saveOBSchedulingInfo(obId, obInfo, obInfoVersion)
        """
        return self.put('/obsBlocks/%d/schedule' % obId, obInfo, etag=version)

# ---------- VISITOR EXECUTION SEQUENCE APIs ----------
    def getExecutionSequence(self, instrument):
        """Retrieve the personal execution sequence for an instrument.

        api: GET /executionSequences/{instrument}

        usage:
           ::

              executionSequence, esVersion = api.getExecutionSequence('UVES')
        """
        return self.get('/executionSequences/%s' % instrument)

    def saveExecutionSequence(self, instrument, executionSequence, version):
        """Change personal visitor execution sequence for an instrument.

        The list should contain objects with an obId in the desired execution
        order and may also contain OBs from delegated observing runs. Only
        visitor mode OBs in status '+', 'A', 'X' or 'S' can be included.

        The obId is the only required field in each list member, but all fields
        sent by a GET request can be included unchanged.

        api: PUT /executionSequences/{instrument}

        usage:
           ::

              executionSequence, esVersion = api.getExecutionSequence('UVES')
              executionSequence, esVersion = api.saveExecutionSequence('UVES',
                 [{ 'obId': ob1 }, { 'obId': ob2 }], esVersion)
              print('OBs in UVES execution sequence', ', '.join(str(e['obId']) for e in executionSequence))
        """
        return self.put('/executionSequences/%s' % instrument, executionSequence, version)

# ---------- INSTRUMENT PACKAGE APIs ----------
    def getInstrumentConstraints(self, instrument, ipVersion):
        """Get all observing constraints that OBs using this instrument package
        must define.

        The information returned is useful for building a user interface and
        also contains allowed values or ranges.

        api: GET /instrumentPackages/{instrument}/{ipVersion}/instrumentConstraints

        usage:
           ::

              insConstraints, _ = api.getInstrumentConstraints('UVES', '98.09')
        """
        return self.get('/instrumentPackages/%s/%s/instrumentConstraints' % (instrument, ipVersion))

    def getTemplateSignatures(self, instrument, ipVersion):
        """List all template signatures that can be attached to OBs using this
        instrument package.

        api: GET /instrumentPackages/{instrument}/{ipVersion}/templateSignatures

        usage:
           ::

              tplSignatures, _ = api.getTemplateSignatures('UVES', '98.09')
        """
        return self.get('/instrumentPackages/%s/%s/templateSignatures' % (instrument, ipVersion))

    def getTemplateSignature(self, instrument, ipVersion, templateName):
        """Get a specific template signature.

        api: GET /instrumentPackages/{instrument}/{ipVersion}/templateSignatures/{templateName}

        usage:
           ::

              tplSignature, _ = api.getTemplateSignature('UVES', '98.08', 'UVES_blue_acq_slit')
        """
        return self.get('/instrumentPackages/%s/%s/templateSignatures/%s' % (instrument, ipVersion, templateName))

    def getInstrumentPackageVersions(self, instrument):
        """Get list of instrument package versions.

        api: GET /instrumentPackages/{instrument}

        usage:
           ::

              ipVersions, _ = api.getInstrumentPackageVersions('UVES')
        """
        return self.get('/instrumentPackages/%s/' % instrument)

# ---------- USER APIs ----------
    def getUser(self):
        """Get information about the logged-in user.

        api: GET /authenticatedUser

        usage:
           ::

              user, _ = api.getUser()
        """
        return self.get('/authenticatedUser')

# ---------- private methods ----------
    def request(self, method, url, data=None, etag=None):
        self.request_count += 1

        # configure request headers
        assert self.access_token is not None
        headers = {
            'Authorization': 'Bearer ' + self.access_token,
            'Accept': 'application/json'
        }
        if data is not None:
            headers['Content-Type'] = 'application/json'
        if etag is not None:
            headers['If-Match'] = etag

        # make request
        url = self.apiUrl + url
        if self.debug:
            print(method, url, data)
        r = requests.request(method, url, headers=headers, data=json.dumps(data))
        content_type = r.headers['Content-Type'].split(';')[0]
        etag = r.headers.get('ETag', None)

        # handle response
        if 200 <= r.status_code < 300:
            if content_type == 'application/json':
                data = r.json()
                return data, etag
            return None, etag
        elif content_type == 'application/json' and 'error' in r.json():
            raise P2Error(r.status_code, method, url, r.json()['error'])
        else:
            raise P2Error(r.status_code, method, url, 'oops unknown error')

    def uploadFile(self, method, url, filename, contentType, etag=None):
        basename = os.path.basename(filename)
        assert self.access_token is not None
        headers = {
            'Authorization': 'Bearer ' + self.access_token,
            'Content-Disposition': 'inline; filename="%s"' % basename,
            'Content-Type': contentType
        }
        if etag is not None:
            headers['If-Match'] = etag
        with open(filename, 'rb') as f:
            r = requests.request(method, self.apiUrl + url, data=f, headers=headers)
            if r.status_code == 201 or r.status_code == 204:
                etag = r.headers.get('ETag', None)
                return basename, etag
            else:
                content_type = r.headers['Content-Type'].split(';')[0]
                if content_type == 'application/json' and 'error' in r.json():
                    raise P2Error(r.status_code, method, url, r.json()['error'])
                else:
                    raise P2Error(r.status_code, method, url, 'oops unknown error')

    def downloadTextFile(self, url, filename):
        assert self.access_token is not None
        headers = {
            'Authorization': 'Bearer ' + self.access_token,
            'Accept': 'text/plain'
        }
        r = requests.request('GET', self.apiUrl + url, stream=True, headers=headers)
        with open(filename, 'wb') as f:
            for line in r.iter_lines():
                f.write(line + b'\n')
        return None, r.headers.get('ETag', None)

    def downloadBinaryFile(self, url, filename):
        assert self.access_token is not None
        headers = {
            'Authorization': 'Bearer ' + self.access_token,
            'Accept': 'text/plain'
        }
        r = requests.request('GET', self.apiUrl + url, stream=True, headers=headers)
        with open(filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=128):
                f.write(chunk)
        return None, None

    def get(self, url):
        return self.request('GET', url)

    def put(self, url, data, etag=None):
        return self.request('PUT', url, data, etag)

    def post(self, url, data=None, etag=None):
        return self.request('POST', url, data, etag)

    def delete(self, url, etag=None):
        return self.request('DELETE', url, etag=etag)
