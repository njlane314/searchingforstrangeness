#! /usr/bin/env python


































































































































































































































































































































































































































































































from __future__ import absolute_import
from __future__ import print_function
import sys, os, stat, subprocess, shutil, json, getpass, uuid, tempfile, hashlib
try:
    import urllib.request as urlrequest
except ImportError:
    import urllib as urlrequest
import larbatch_posix
import threading
try:
    import queue
except ImportError:
    import Queue as queue
from xml.dom.minidom import parse
import project_utilities
from project_modules.projectdef import ProjectDef
from project_modules.projectstatus import ProjectStatus
from project_modules.batchstatus import BatchStatus
from project_modules.jobsuberror import JobsubError
from project_modules.ifdherror import IFDHError
import larbatch_utilities
from larbatch_utilities import convert_str
from larbatch_utilities import convert_bytes
import samweb_cli

samweb = None
extractor_dict = None





def import_samweb():



    global samweb
    global extractor_dict
    global expMetaData


    if samweb == None:
        samweb = project_utilities.samweb()
        from extractor_dict import expMetaData



def docleanx(projects, projectname, stagename, clean_descendants = True):
    print(projectname, stagename)











    uid = os.getuid()
    euid = os.geteuid()
    cleaned_bookdirs = []



    done_cleaning = False
    while not done_cleaning:

        cleaned_something = False



        for project in projects:
            for stage in project.stages:

                clean_this_stage = False



                if not stage.bookdir in cleaned_bookdirs:



                    if (projectname == '' or project.name == projectname) and \
                       (stagename == '' or stage.name == stagename):

                        clean_this_stage = True




                    elif clean_descendants and stage.inputlist != '' and \
                         os.path.dirname(stage.inputlist) in cleaned_bookdirs:

                        clean_this_stage = True



                    if clean_this_stage:
                        cleaned_something = True
                        cleaned_bookdirs.append(stage.bookdir)

                        print('Clean project %s, stage %s' % (project.name, stage.name))



                        if larbatch_posix.exists(stage.outdir):
                            dir_uid = larbatch_posix.stat(stage.outdir).st_uid
                            if dir_uid == uid or dir_uid == euid:
                                print('Clean directory %s.' % stage.outdir)
                                larbatch_posix.rmtree(stage.outdir)
                            else:
                                raise RuntimeError('Owner mismatch, delete %s manually.' % stage.outdir)



                        if larbatch_posix.exists(stage.logdir):
                            dir_uid = larbatch_posix.stat(stage.logdir).st_uid
                            if dir_uid == uid or dir_uid == euid:
                                print('Clean directory %s.' % stage.logdir)
                                larbatch_posix.rmtree(stage.logdir)
                            else:
                                raise RuntimeError('Owner mismatch, delete %s manually.' % stage.logdir)



                        if larbatch_posix.exists(stage.workdir):
                            dir_uid = larbatch_posix.stat(stage.workdir).st_uid
                            if dir_uid == uid or dir_uid == euid:
                                print('Clean directory %s.' % stage.workdir)
                                larbatch_posix.rmtree(stage.workdir)
                            else:
                                raise RuntimeError('Owner mismatch, delete %s manually.' % stage.workdir)



                        if larbatch_posix.exists(stage.bookdir):
                            dir_uid = larbatch_posix.stat(stage.bookdir).st_uid
                            if dir_uid == uid or dir_uid == euid:
                                print('Clean directory %s.' % stage.bookdir)
                                larbatch_posix.rmtree(stage.bookdir)
                            else:
                                raise RuntimeError('Owner mismatch, delete %s manually.' % stage.bookdir)

        done_cleaning = not cleaned_something



    return



def dostatus(projects):



    project_utilities.test_token()




    prjs = projects
    if type(projects) != type([]) and type(projects) != type(()):
        prjs = [projects]

    project_status = ProjectStatus(prjs)
    batch_status = BatchStatus(prjs)

    for project in prjs:

        print('\nProject %s:' % project.name)



        for stage in project.stages:

            stagename = stage.name
            stage_status = project_status.get_stage_status(stagename)
            b_stage_status = batch_status.get_stage_status(stagename)
            if stage_status.exists:
                print('\nStage %s: %d art files, %d events, %d analysis files, %d errors, %d missing files.' % (
                    stagename, stage_status.nfile, stage_status.nev, stage_status.nana,
                    stage_status.nerror, stage_status.nmiss))
            else:
                print('\nStage %s output directory does not exist.' % stagename)
            print('Stage %s batch jobs: %d idle, %d running, %d held, %d other.' % (
                stagename, b_stage_status[0], b_stage_status[1], b_stage_status[2], b_stage_status[3]))
    return




def find_projects(element, check=True):

    projects = []




    if element.nodeName == 'project':
        default_input_by_stage = {}
        project = ProjectDef(element, '', default_input_by_stage, check=check)
        projects.append(project)

    else:




        default_input = ''
        default_input_by_stage = {}
        subelements = element.getElementsByTagName('project')
        for subelement in subelements:
            project = ProjectDef(subelement, default_input, default_input_by_stage, check=check)
            projects.append(project)
            for stage in project.stages:
                stage_list = os.path.join(stage.bookdir, 'files.list')
                default_input_by_stage[stage.name] = stage_list
                default_input = stage_list



    return projects




def get_projects(xmlfile, check=True):



    if xmlfile in get_projects.cache:
        return get_projects.cache[xmlfile]



    if xmlfile == '-':
        xml = sys.stdin
    elif xmlfile.find(':') < 0:
        xml = open(xmlfile)
    else:
        xml = urlrequest.urlopen(xmlfile)
    doc = parse(xml)



    root = doc.documentElement



    projects = find_projects(root, check=check)



    get_projects.cache[xmlfile] = projects



    return projects



get_projects.cache = {}




def select_project(projects, projectname, stagename):

    for project in projects:
        if projectname == '' or projectname == project.name:
            for stage in project.stages:
                if stagename == '' or stagename == stage.name:
                    return project



    return None




def get_project(xmlfile, projectname='', stagename='', check=True):
    projects = get_projects(xmlfile, check=check)
    project = select_project(projects, projectname, stagename)
    return project



def next_stage(projects, stagename, circular=False):



    found = False
    for project in projects:



        for stage in project.stages:
            if found:
                return stage
            if stage.name == stagename:
                found = True



    if circular and len(projects) > 0 and len(projects[0].stages) > 0:
        return projects[0].stages[0]



    return None



def previous_stage(projects, stagename, circular=False):



    result = None
    if circular and len(projects) > 0 and len(projects[-1].stages) > 0:
        result = projects[-1].stages[-1]



    for project in projects:



        for stage in project.stages:
            if stage.name == stagename:
                return result
            result = stage



    return result




def get_pubs_stage(xmlfile, projectname, stagename, run, subruns, version=None):
    projects = get_projects(xmlfile)
    project = select_project(projects, projectname, stagename)
    if project == None:
        raise RuntimeError('No project selected for projectname=%s, stagename=%s' % (
            projectname, stagename))
    stage = project.get_stage(stagename)
    if stage == None:
        raise RuntimeError('No stage selected for projectname=%s, stagename=%s' % (
            projectname, stagename))
    get_projects.cache = {}
    stage.pubsify_input(run, subruns, version)
    stage.pubsify_output(run, subruns, version)
    get_projects.cache = {}
    return project, stage









def check_root_file(path, logdir):

    result = (-2, '')
    json_ok = False
    md = []



    if not larbatch_posix.exists(path):
        return result



    json_path = os.path.join(logdir, os.path.basename(path) + '.json')
    if larbatch_posix.exists(json_path):



        try:
            lines = larbatch_posix.readlines(json_path)
            s = ''
            for line in lines:
                s = s + line



            md = json.loads(s)



            result = (-1, '')



            if len(list(md.keys())) > 0:
                nevroot = -1
                stream = ''
                if 'events' in md:
                    nevroot = int(md['events'])
                if 'data_stream' in md:
                    stream = md['data_stream']
                result = (nevroot, stream)
            json_ok = True
        except:
            result = (-2, '')
    return result




def check_root(outdir, logdir, data_file_types):














    nev = -1
    roots = []
    hists = []

    print('Checking root files in directory %s.' % outdir)
    filenames = larbatch_posix.listdir(outdir)
    for filename in filenames:
        name, ext = os.path.splitext(filename)
        if len(ext) > 0 and ext[1:] in data_file_types:
            path = os.path.join(outdir, filename)
            nevroot, stream = check_root_file(path, logdir)
            if nevroot >= 0:
                if nev < 0:
                    nev = 0
                nev = nev + nevroot
                roots.append((os.path.join(outdir, filename), nevroot, stream))

            elif nevroot == -1:



                hists.append(os.path.join(outdir, filename))

            else:




                print('Warning: File %s in directory %s is not a valid root file.' % (filename, outdir))



    return (nev, roots, hists)




def get_input_files(stage):





    result = []
    if stage.inputfile != '':
        result.append(stage.inputfile)

    elif stage.inputlist != '' and larbatch_posix.exists(stage.inputlist):
        try:
            input_filenames = larbatch_posix.readlines(stage.inputlist)
            for line in input_filenames:
                words = line.split()
                result.append(words[0])
        except:
            pass

    elif stage.inputdef != '':
        import_samweb()
        result = samweb.listFiles(defname=stage.inputdef)



    return result



def doshorten(stage):



    untarlog(stage)



    for out_subpath, subdirs, files in larbatch_posix.walk(stage.outdir):



        if len(subdirs) != 0:
            continue

        subdir = os.path.relpath(out_subpath, stage.outdir)
        log_subpath = os.path.join(stage.bookdir, subdir)

        for file in files:
            if file[-5:] == '.root':
                if len(file) >= 200:



                    file_path = os.path.join(out_subpath, file)
                    shortfile = file[:150] + str(uuid.uuid4()) + '.root'
                    shortfile_path = os.path.join(out_subpath, shortfile)
                    print('%s\n->%s\n' % (file_path, shortfile_path))
                    larbatch_posix.rename(file_path, shortfile_path)



                    json_path = os.path.join(log_subpath, file + '.json')
                    if larbatch_posix.exists(json_path):
                        shortjson = shortfile + '.json'
                        shortjson_path = os.path.join(log_subpath, shortjson)
                        print('%s\n->%s\n' % (json_path, shortjson_path))
                        larbatch_posix.rename(json_path, shortjson_path)

    return



def untarlog(stage):



    for log_subpath, subdirs, files in larbatch_posix.walk(stage.logdir):



        if len(subdirs) != 0:
            continue
        subdir = os.path.relpath(log_subpath, stage.logdir)
        if subdir == '.':
            continue
        book_subpath = os.path.join(stage.bookdir, subdir)
        for file in files:
            if file.startswith('log') and file.endswith('.tar'):
                src = '%s/%s' % (log_subpath, file)
                dst = '%s/%s' % (book_subpath, file)
                flag = '%s.done' % dst



                if dst != src and not larbatch_posix.exists(flag):



                    print('Copying tarball %s into %s' % (src, book_subpath))
                    if not larbatch_posix.isdir(book_subpath):
                        larbatch_posix.makedirs(book_subpath)
                    larbatch_posix.copy(src, dst)



                if not larbatch_posix.exists(flag):



                    print('Extracting tarball %s' % dst)
                    jobinfo = subprocess.Popen(['tar','-xf', dst, '-C', book_subpath,
                                                '--exclude=beam*.dat',
                                                '--exclude=beam*.info',
                                                '--exclude=core*',
                                                '--exclude=*.db',
                                                '--exclude=*.sh',
                                                '--exclude=*.py*',
                                                '--exclude=*.tar'],
                                               stdout=subprocess.PIPE,
                                               stderr=subprocess.PIPE)
                    jobout, joberr = jobinfo.communicate()
                    jobout = convert_str(jobout)
                    joberr = convert_str(joberr)
                    rc = jobinfo.poll()
                    if rc != 0:
                        print(jobout)
                        print(joberr)
                        print('Failed to extract log tarball in %s' % dst)

                    else:



                        f = larbatch_posix.open(flag, 'w')
                        f.write('\n')
                        f.close()



                        if dst != src:
                            larbatch_posix.remove(dst)

    return



def docheck(project, stage, ana, quick=False):


















































    untarlog(stage)



    if quick == 1 and not ana:
        return doquickcheck(project, stage, ana)

    stage.checkinput()



    if not larbatch_posix.exists(stage.outdir):
        print('Output directory %s does not exist.' % stage.outdir)
        return 1
    if not larbatch_posix.exists(stage.bookdir):
        print('Log directory %s does not exist.' % stage.bookdir)
        return 1

    import_samweb()
    has_metadata = project.file_type != '' or project.run_type != ''
    has_input = stage.inputfile != '' or stage.inputlist != '' or stage.inputdef != ''
    print('Checking directory %s' % stage.bookdir)



    nev_tot = 0
    nroot_tot = 0



    procmap = {}
    processes = []
    filesana = []
    sam_projects = []
    cpids = []
    uris = []
    bad_workers = []


    for log_subpath, subdirs, files in larbatch_posix.walk(stage.bookdir):



        if len(subdirs) != 0:
            continue

        subdir = os.path.relpath(log_subpath, stage.bookdir)
        if subdir == '.':
            continue
        out_subpath = os.path.join(stage.outdir, subdir)
        dirok = project_utilities.fast_isdir(log_subpath)



        if dirok and log_subpath[-6:] == '_start':
            filename = os.path.join(log_subpath, 'sam_project.txt')
            if larbatch_posix.exists(filename):
                sam_project = larbatch_posix.readlines(filename)[0].strip()
                if sam_project != '' and not sam_project in sam_projects:
                    sam_projects.append(sam_project)



        if dirok and not subdir[-6:] == '_start' and not subdir[-5:] == '_stop' \
                and not subdir == 'log':

            bad = 0



            if not project_utilities.fast_isdir(out_subpath):
                print('No output directory corresponding to subdirectory %s.' % subdir)
                bad = 1



            if not bad:
                stat_filename = os.path.join(log_subpath, 'lar.stat')
                if larbatch_posix.exists(stat_filename):
                    status = 0
                    try:
                        status = int(larbatch_posix.readlines(stat_filename)[0].strip())
                        if status != 0:
                            print('Job in subdirectory %s ended with non-zero exit status %d.' % (
                                subdir, status))
                            bad = 1
                    except:
                        print('Bad file lar.stat in subdirectory %s.' % subdir)
                        bad = 1



            if not bad:
                nev = 0
                roots = []
                nev, roots, subhists = check_root(out_subpath, log_subpath, stage.datafiletypes)
                if not ana:
                    if len(roots) == 0 or nev < 0:
                        print('Problem with root file(s) in subdirectory %s.' % subdir)
                        bad = 1
                elif nev < -1 or len(subhists) == 0:
                    print('Problem with analysis root file(s) in subdirectory %s.' % subdir)
                    bad = 1




            if not bad and has_metadata:
                for root in roots:
                    rootname = os.path.basename(root[0])
                    for s in list(procmap.keys()):
                        oldroots = procmap[s]
                        for oldroot in oldroots:
                            oldrootname = os.path.basename(oldroot[0])
                            if rootname == oldrootname:
                                print('Duplicate filename %s in subdirectory %s' % (rootname,
                                                                                    subdir))
                                olddir = os.path.basename(os.path.dirname(oldroot[0]))
                                print('Previous subdirectory %s' % olddir)
                                bad = 1



            if not bad and has_metadata:
                for root in roots:
                    rootname = os.path.basename(root[0])
                    if len(rootname) >= 200:
                        print('Filename %s in subdirectory %s is longer than 200 characters.' % (
                            rootname, subdir))
                        bad = 1




            if not bad and stage.inputdef != '':
                filename1 = os.path.join(log_subpath, 'sam_project.txt')
                if not larbatch_posix.exists(filename1):
                    print('Could not find file sam_project.txt')
                    bad = 1
                filename2 = os.path.join(log_subpath, 'cpid.txt')
                if not larbatch_posix.exists(filename2):
                    print('Could not find file cpid.txt')
                    bad = 1
                if not bad:
                    sam_project = larbatch_posix.readlines(filename1)[0].strip()
                    if not sam_project in sam_projects:
                        sam_projects.append(sam_project)
                    cpid = larbatch_posix.readlines(filename2)[0].strip()
                    if not cpid in cpids:
                        cpids.append(cpid)




            if not bad and (stage.inputlist !='' or stage.inputfile != ''):
                filename = os.path.join(log_subpath, 'transferred_uris.list')
                if not larbatch_posix.exists(filename):
                    print('Could not find file transferred_uris.list')
                    bad = 1
                if not bad:
                    lines = larbatch_posix.readlines(filename)
                    for line in lines:
                        uri = line.strip()
                        if uri != '':
                            uris.append(uri)




            if not has_input:
                subdir_split = subdir.split('_')
                if len(subdir_split) > 1:
                    process = int(subdir_split[1])
                    if process in processes:
                        print('Duplicate process number')
                        bad = 1
                    else:
                        processes.append(process)



            if not bad:
                procmap[subdir] = roots



                filesana.extend(subhists)



                nev_tot = nev_tot + nev
                nroot_tot = nroot_tot + len(roots)



            if bad:
                bad_workers.append(subdir)



            if bad:
                print('Bad subdirectory %s.' % subdir)










    contents = larbatch_posix.listdir(stage.bookdir)
    if len(contents) == 0:
        print('Directory %s may be dead.' % stage.bookdir)
        print('Returning error status without creating any bookkeeping files.')
        return 1



    filelistname = os.path.join(stage.bookdir, 'files.list')
    filelist = safeopen(filelistname)

    eventslistname = os.path.join(stage.bookdir, 'events.list')
    eventslist = safeopen(eventslistname)

    badfilename = os.path.join(stage.bookdir, 'bad.list')
    badfile = safeopen(badfilename)

    missingfilesname = os.path.join(stage.bookdir, 'missing_files.list')
    missingfiles = safeopen(missingfilesname)

    filesanalistname = os.path.join(stage.bookdir, 'filesana.list')
    filesanalist = safeopen(filesanalistname)

    urislistname = os.path.join(stage.bookdir, 'transferred_uris.list')
    urislist = safeopen(urislistname)




    nproc = 0
    streams = {}
    nfile = 0
    for s in list(procmap.keys()):
        nproc = nproc + 1
        for root in procmap[s]:
            nfile = nfile + 1
            filelist.write('%s\n' % root[0])
            eventslist.write('%s %d\n' % root[:2])
            stream = root[2]
            if stream != '':
                if stream not in streams:
                    streamlistname = os.path.join(stage.bookdir, 'files_%s.list' % stream)
                    streams[stream] = safeopen(streamlistname)
                streams[stream].write('%s\n' % root[0])



    nerror = 0
    for bad_worker in bad_workers:
        badfile.write('%s\n' % bad_worker)
        nerror = nerror + 1



    nmiss = 0
    if stage.inputdef == '' and not stage.pubs_output:
        input_files = get_input_files(stage)
        if len(input_files) > 0:
            missing_files = list(set(input_files) - set(uris))
            for missing_file in missing_files:
                missingfiles.write('%s\n' % missing_file)
                nmiss = nmiss + 1
        else:
            nmiss = stage.num_jobs - len(procmap)
            for n in range(nmiss):
                missingfiles.write('/dev/null\n')




    for hist in filesana:
        filesanalist.write('%s\n' % hist)



    for uri in uris:
        urislist.write('%s\n' % uri)



    if ana:
        print("%d processes completed successfully." % nproc)
        print("%d total good histogram files." % len(filesana))
    else:
        print("%d total good events." % nev_tot)
        print("%d total good root files." % nroot_tot)
        print("%d total good histogram files." % len(filesana))



    filelist.close()
    if nfile == 0:
        project_utilities.addLayerTwo(filelistname)
    eventslist.close()
    if nfile == 0:
        project_utilities.addLayerTwo(eventslistname)
    if nerror == 0:
        badfile.write('\n')
    badfile.close()
    if nmiss == 0:
        missingfiles.write('\n')
    missingfiles.close()
    filesanalist.close()
    if len(filesana) == 0:
        project_utilities.addLayerTwo(filesanalistname)
    if len(uris) == 0:
        urislist.write('\n')
    urislist.close()
    for stream in list(streams.keys()):
        streams[stream].close()



    if stage.inputdef != '' and not stage.pubs_input:



        sam_projects_filename = os.path.join(stage.bookdir, 'sam_projects.list')
        sam_projects_file = safeopen(sam_projects_filename)
        for sam_project in sam_projects:
            sam_projects_file.write('%s\n' % sam_project)
        sam_projects_file.close()
        if len(sam_projects) == 0:
            project_utilities.addLayerTwo(sam_projects_filename)



        cpids_filename = os.path.join(stage.bookdir, 'cpids.list')
        cpids_file = safeopen(cpids_filename)
        for cpid in cpids:
            cpids_file.write('%s\n' % cpid)
        cpids_file.close()
        if len(cpids) == 0:
            project_utilities.addLayerTwo(cpids_filename)



        cpids_list = ''
        sep = ''
        for cpid in cpids:
            cpids_list = cpids_list + '%s%s' % (sep, cpid)
            sep = ','
        if cpids_list != '':
            dim = 'consumer_process_id %s and consumed_status consumed' % cpids_list
            import_samweb()
            nconsumed = samweb.countFiles(dimensions=dim)
        else:
            nconsumed = 0



        if cpids_list != '':
            udim = '(defname: %s) minus (%s)' % (stage.inputdef, dim)
        else:
            udim = 'defname: %s' % stage.inputdef
        nunconsumed = samweb.countFiles(dimensions=udim)
        nerror = nerror + nunconsumed



        print('%d sam projects.' % len(sam_projects))
        print('%d successful consumer process ids.' % len(cpids))
        print('%d files consumed.' % nconsumed)
        print('%d files not consumed.' % nunconsumed)



        for sam_project in sam_projects:
            print('\nChecking sam project %s' % sam_project)
            import_samweb()
            url = samweb.findProject(sam_project, project_utilities.get_experiment())
            if url != '':
                result = samweb.projectSummary(url)
                nd = 0
                nc = 0
                nf = 0
                nproc = 0
                nact = 0
                if 'processes' in result:
                    processes = result['processes']
                    for process in processes:
                        nproc = nproc + 1
                        if 'status' in process:
                            if process['status'] == 'active':
                                nact = nact + 1
                        if 'counts' in process:
                            counts = process['counts']
                            if 'delivered' in counts:
                                nd = nd + counts['delivered']
                            if 'consumed' in counts:
                                nc = nc + counts['consumed']
                            if 'failed' in counts:
                                nf = nf + counts['failed']
                print('Status: %s' % result['project_status'])
                print('%d total processes' % nproc)
                print('%d active processes' % nact)
                print('%d files in snapshot' % result['files_in_snapshot'])
                print('%d files delivered' % (nd + nc))
                print('%d files consumed' % nc)
                print('%d files failed' % nf)
                print()



    checkfilename = os.path.join(stage.bookdir, 'checked')
    checkfile = safeopen(checkfilename)
    checkfile.write('\n')
    checkfile.close()
    project_utilities.addLayerTwo(checkfilename)

    if stage.inputdef == '' or stage.pubs_input:
        print('%d processes with errors.' % nerror)
        print('%d missing files.' % nmiss)
    else:
        print('%d unconsumed files.' % nerror)




    result = 0
    if nerror != 0:
        result = 1
    if not ana and nroot_tot == 0:
        result = 1
    if len(procmap) == 0:
        result = 1
    return result

def doquickcheck(project, stage, ana):


    if not larbatch_posix.isdir(stage.outdir):
        print('Output directory %s does not exist.' % stage.outdir)
        return 1

    if not larbatch_posix.isdir(stage.bookdir):
        print('Log directory %s does not exist.' % stage.bookdir)
        return 1

    print('Checking directory %s' % stage.bookdir)



    goodFiles        = []
    goodAnaFiles     = []
    eventLists       = []
    badLists         = []
    anaFiles         = []
    transferredFiles = []
    streamLists      = {}

    sam_projects     = []
    cpids            = []

    goodLogDirs      = set()
    nErrors = 0

    for log_subpath, subdirs, files in larbatch_posix.walk(stage.bookdir):



        if len(subdirs) != 0:
            continue


        if log_subpath[-6:] == '_start' or log_subpath[-5:] == '_stop':
            filename = os.path.join(log_subpath, 'sam_project.txt')
            if larbatch_posix.exists(filename):
                sam_project = larbatch_posix.readlines(filename)[0].strip()
                if sam_project != '' and not sam_project in sam_projects:
                    sam_projects.append(sam_project)
            continue


        print('Doing quick check of directory %s.' % log_subpath)

        subdir = os.path.relpath(log_subpath, stage.bookdir)

        out_subpath = os.path.join(stage.outdir, subdir)
        dirok = project_utilities.fast_isdir(log_subpath)




        validateOK = 1

        missingfilesname = os.path.join(log_subpath, 'missing_files.list')



        try:

            missingfiles = project_utilities.saferead(missingfilesname)

        except:
            print('Cannot open file: %s' % missingfilesname)
            validateOK = 0


        if validateOK == 1 and len(missingfiles) == 0:
            print('%s exists, but is empty' % missingfilesname)
            validateOK = 0


        if validateOK == 1:
            line = missingfiles[0]
            line = line.strip('\n')
            if( int(line) != 0 ):
                validateOK = 0



        if validateOK != 1:
            nErrors += 1
            continue







        if stage.inputdef != '':

            filename1 = os.path.join(log_subpath, 'sam_project.txt')
            if not larbatch_posix.exists(filename1):
                print('Could not find file sam_project.txt')
                nErrors += 1
            else:
                sam_project = larbatch_posix.readlines(filename1)[0].strip()
                if not sam_project in sam_projects:
                    sam_projects.append(sam_project)

            filename2 = os.path.join(log_subpath, 'cpid.txt')
            if not larbatch_posix.exists(filename2):
                print('Could not find file cpid.txt')
                nErrors += 1
            else:
                cpid = larbatch_posix.readlines(filename2)[0].strip()
                if not cpid in cpids:
                    cpids.append(cpid)

        filelistsrc = os.path.join(log_subpath, 'files.list')
        tmpArray = scan_file(filelistsrc)

        if( tmpArray == [ -1 ] ):
            nErrors += 1
        else:
            goodFiles.extend(tmpArray)

        fileanalistsrc = os.path.join(log_subpath, 'filesana.list')
        tmpArray = scan_file(fileanalistsrc)

        if( not tmpArray == [ -1 ] ):
            goodAnaFiles.extend(tmpArray)

        eventlistsrc = os.path.join(log_subpath, 'events.list')

        tmpArray = scan_file(eventlistsrc)

        if( tmpArray == [ -1 ] ):
            nErrors += 1
        else:
            eventLists.extend(tmpArray)


        badfilesrc = os.path.join(log_subpath, 'bad.list')


        tmpArray = scan_file(badfilesrc)


        if( tmpArray == [ -1 ] ):
            pass
        else:
            badLists.extend(tmpArray)

        '''
        missingfilesrc  = os.path.join(log_subpath, 'missing_files.list')

        tmpArray = scan_file(missingfilesrc)

        if( tmpArray == [ -1 ] ):
            nErrors += 1
        else:
            missingLists.extend(tmpArray)
        '''











        urislistsrc = os.path.join(log_subpath, 'transferred_uris.list')

        tmpArray = scan_file(urislistsrc)


        if( tmpArray == [ -1 ] ):
            pass
        else:
            transferredFiles.extend(tmpArray)

        streamList = larbatch_posix.listdir(log_subpath)

        for stream in streamList:
            if( stream[:6] != "files_" ):
                continue
            streamfilesrc = os.path.join(log_subpath, stream)

            tmpArray = scan_file(streamfilesrc)
            if( tmpArray == [ -1 ] ):
                nErrors += 1
            else:
                if(streamLists.get(stream, "empty") == "empty" ):
                    streamLists[stream] = tmpArray
                else:
                    streamLists[stream].extend(tmpArray)

        if validateOK == 1:
            goodLogDirs.add(log_subpath)

    checkfilename = os.path.join(stage.bookdir, 'checked')
    checkfile = safeopen(checkfilename)
    checkfile.write('\n')
    checkfile.close()


    filelistdest = os.path.join(stage.bookdir, 'files.list')
    if larbatch_posix.exists(filelistdest):

        larbatch_posix.remove(filelistdest)
    if len(goodLogDirs) == 1:
        src = '%s/files.list' % goodLogDirs.copy().pop()

        larbatch_posix.symlink(src, filelistdest)
    else:

        inputList = safeopen(filelistdest)
        for goodFile in goodFiles:

            inputList.write("%s\n" % goodFile)
        inputList.close()
    if len(goodFiles) == 0:
        project_utilities.addLayerTwo(filelistdest)


    fileanalistdest = os.path.join(stage.bookdir, 'filesana.list')
    if larbatch_posix.exists(fileanalistdest):

        larbatch_posix.remove(fileanalistdest)
    if len(goodLogDirs) == 1:
        src = '%s/filesana.list' % goodLogDirs.copy().pop()

        larbatch_posix.symlink(src, fileanalistdest)
    else:

        anaList = safeopen(fileanalistdest)
        for goodAnaFile in goodAnaFiles:

            anaList.write("%s\n" % goodAnaFile)
        anaList.close()
        if len(goodAnaFiles) == 0:
            project_utilities.addLayerTwo(fileanalistdest)


    eventlistdest = os.path.join(stage.bookdir, 'events.list')
    if larbatch_posix.exists(eventlistdest):

        larbatch_posix.remove(eventlistdest)
    if len(goodLogDirs) == 1:
        src = '%s/events.list' % goodLogDirs.copy().pop()

        larbatch_posix.symlink(src, eventlistdest)
    else:

        eventsOutList = safeopen(eventlistdest)
        for event in eventLists:

            eventsOutList.write("%s\n" % event)
        eventsOutList.close()
        if len(eventLists) == 0:
            project_utilities.addLayerTwo(eventlistdest)


    if(len(badLists) > 0):
        badlistdest = os.path.join(stage.bookdir, 'bad.list')
        badOutList = safeopen(badlistdest)
        for bad in badLists:
            badOutList.write("%s\n" % bad)
        badOutList.close()



    missing_files = []
    if stage.inputdef == '' and not stage.pubs_output:
        input_files = get_input_files(stage)
        if len(input_files) > 0:
            missing_files = list(set(input_files) - set(transferredFiles))

    if len(missing_files) > 0:
        missinglistdest = os.path.join(stage.bookdir, 'missing_files.list')
        missingOutList = safeopen(missinglistdest)
        for missing in missing_files:
            missingOutList.write("%s\n" % missing)
        missingOutList.close()



    urilistdest = os.path.join(stage.bookdir, 'transferred_uris.list')
    if larbatch_posix.exists(urilistdest):

        larbatch_posix.remove(urilistdest)
    if len(goodLogDirs) == 1 and len(transferredFiles) > 0:
        src = '%s/transferred_uris.list' % goodLogDirs.copy().pop()

        larbatch_posix.symlink(src, urilistdest)
    else:

        uriOutList  = safeopen(urilistdest)
        for uri in transferredFiles:

            uriOutList.write("%s\n" % uri)
        uriOutList.close()
        if len(transferredFiles) == 0:
            project_utilities.addLayerTwo(urilistdest)

    if stage.inputdef != '':
        samprojectdest = os.path.join(stage.bookdir, 'sam_projects.list')
        if larbatch_posix.exists(samprojectdest):

            larbatch_posix.remove(samprojectdest)
        if len(goodLogDirs) == 1:
            src = '%s/sam_project.txt' % goodLogDirs.copy().pop()

            larbatch_posix.symlink(src, samprojectdest)
        else:

            samprojectfile = safeopen(samprojectdest)
            for sam in sam_projects:
                samprojectfile.write("%s\n" % sam)
            samprojectfile.close()
            if len(sam_projects) == 0:
                project_utilities.addLayerTwo(samprojectdest)

        cpiddest = os.path.join(stage.bookdir, 'cpids.list')
        if larbatch_posix.exists(cpiddest):

            larbatch_posix.remove(cpiddest)
        if len(goodLogDirs) == 1:
            src = '%s/cpid.txt' % goodLogDirs.copy().pop()

            larbatch_posix.symlink(src, cpiddest)
        else:

            cpidfile = safeopen(cpiddest)
            for cp in cpids:
                cpidfile.write("%s \n" % cp)
            cpidfile.close()
            if len(cpids) == 0:
                project_utilities.addLayerTwo(cpiddest)


    for stream in streamLists:
        streamdest = os.path.join(stage.bookdir, stream)
        if larbatch_posix.exists(streamdest):

            larbatch_posix.remove(streamdest)
        if len(goodLogDirs) == 1:
            src = '%s/%s' % (goodLogDirs.copy().pop(), stream)

            larbatch_posix.symlink(src, streamdest)
        else:

            streamOutList = safeopen(streamdest)
            for line in streamLists[stream]:
                streamOutList.write("%s\n" % line)
            streamOutList.close()
            if len(streamLists[stream]) == 0:
                project_utilities.addLayerTwo(streamdest)





    print('Number of errors = %d' % nErrors)

    return nErrors



def dofetchlog(project, stage):











    stage.checkinput()
    stage.checkdirs()




    logids = []
    for dirpath, dirnames, filenames in larbatch_posix.walk(stage.bookdir):
        for filename in filenames:
            if filename == 'env.txt':









                logid = ''
                envpath = os.path.join(dirpath, filename)
                vars = larbatch_posix.readlines(envpath)



                for var in vars:
                    varsplit = var.split('=', 1)
                    name = varsplit[0].strip()
                    if name == 'JOBSUBPARENTJOBID':
                        logid = varsplit[1].strip()




                        logsplit = logid.split('@', 1)
                        cluster_process = logsplit[0]
                        server = logsplit[1]
                        cluster = cluster_process.split('.', 1)[0]
                        logid = cluster + '.0' + '@' + server
                        logids.append(logid)
                        break



                if logid == '':
                    for var in vars:
                        varsplit = var.split('=', 1)
                        name = varsplit[0].strip()
                        if name == 'JOBSUBJOBID':
                            logid = varsplit[1].strip()




                            logsplit = logid.split('@', 1)
                            cluster_process = logsplit[0]
                            server = logsplit[1]
                            cluster = cluster_process.split('.', 1)[0]
                            logid = cluster + '.0' + '@' + server
                            logids.append(logid)
                            break



    if len(logids) > 0:



        logdir = os.path.join(stage.bookdir, 'log')
        if larbatch_posix.exists(logdir):
            larbatch_posix.rmtree(logdir)
        larbatch_posix.mkdir(logdir)



        for logid in set(logids):





            print('Fetching log files for id %s' % logid)
            command = ['jobsub_fetchlog']
            if project.server != '-' and project.server != '':
                command.append('--jobsub-server=%s' % project.server)
            command.append('--jobid=%s' % logid)
            command.append('--dest-dir=%s' % logdir)
            jobinfo = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            jobout, joberr = jobinfo.communicate()
            jobout = convert_str(jobout)
            joberr = convert_str(joberr)
            rc = jobinfo.poll()
            if rc != 0:
                raise JobsubError(command, rc, jobout, joberr)

        return 0

    else:







        print('Failed to fetch log files.')
        return 1






def docheck_declarations(logdir, outdir, declare, ana=False):



    result = 0



    import_samweb()



    roots = []
    listname = 'files.list'
    if ana:
        listname = 'filesana.list'
    fnlist = os.path.join(logdir, listname)
    if larbatch_posix.exists(fnlist):
        roots = larbatch_posix.readlines(fnlist)
    else:
        raise RuntimeError('No %s file found %s, run project.py --check' % (listname, fnlist))

    for root in roots:
        path = root.strip()
        fn = os.path.basename(path)
        dirpath = os.path.dirname(path)
        dirname = os.path.relpath(dirpath, outdir)



        has_metadata = False
        try:
            md = samweb.getMetadata(filenameorid=fn)
            has_metadata = True
        except samweb_cli.exceptions.FileNotFound:
            pass



        if has_metadata:
            print('Metadata OK: %s' % fn)
        else:
            if declare:
                print('Declaring: %s' % fn)
                jsonfile = os.path.join(logdir, os.path.join(dirname, fn)) + '.json'
                mdjson = {}
                if larbatch_posix.exists(jsonfile):
                    mdlines = larbatch_posix.readlines(jsonfile)
                    mdtext = ''
                    for line in mdlines:
                        mdtext = mdtext + line
                    try:
                        md = json.loads(mdtext)
                        mdjson = md
                    except:
                        pass
                md = {}
                if ana:
                    md = mdjson
                else:
                    expSpecificMetaData = expMetaData(os.environ['SAM_EXPERIMENT'],larbatch_posix.root_stream(path))
                    md = expSpecificMetaData.getmetadata(mdjson)
                if len(md) > 0:
                    project_utilities.test_token()




                    try:
                        samweb.declareFile(md=md)
                    except:



                        print('SAM declare failed.')
                        result = 1

                else:
                    print('No sam metadata found for %s.' % fn)
            else:
                print('Not declared: %s' % fn)
                result = 1

    return result



def dotest_declarations(dim):



    import_samweb()



    result = samweb.listFilesSummary(dimensions=dim)
    for key in list(result.keys()):
        print('%s: %s' % (key, result[key]))

    return 0





def docheck_definition(defname, dim, define):



    result = 0



    if defname == '':
        return result



    import_samweb()



    def_exists = False
    try:
        desc = samweb.descDefinition(defname=defname)
        def_exists = True
    except samweb_cli.exceptions.DefinitionNotFound:
        pass



    if def_exists:
        print('Definition already exists: %s' % defname)
    else:
        if define:
            print('Creating definition %s.' % defname)
            project_utilities.test_token()
            samweb.createDefinition(defname=defname, dims=dim)
        else:
            result = 1
            print('Definition should be created: %s' % defname)

    return result



def dotest_definition(defname):



    import_samweb()



    result = samweb.listFilesSummary(defname=defname)
    for key in list(result.keys()):
        print('%s: %s' % (key, result[key]))

    return 0



def doundefine(defname):

    if defname == '':
        return 1



    import_samweb()



    def_exists = False
    try:
        desc = samweb.descDefinition(defname=defname)
        def_exists = True
    except samweb_cli.exceptions.DefinitionNotFound:
        pass



    if def_exists:
        print('Deleting definition: %s' % defname)
        project_utilities.test_token()
        samweb.deleteDefinition(defname=defname)
    else:
        print('No such definition: %s' % defname)

    return 0




def docheck_locations(dim, outdir, add, clean, remove, upload):

    if add:
        print('Adding disk locations.')
    elif clean:
        print('Cleaning disk locations.')
    elif remove:
        print('Removing disk locations.')
    elif upload:
        print('Uploading to FTS.')
    else:
        print('Checking disk locations.')



    import_samweb()



    filelist = samweb.listFiles(dimensions=dim, stream=False)



    disk_dict = {}
    for filename in filelist:
        disk_dict[filename] = []
    for out_subpath, subdirs, files in larbatch_posix.walk(outdir):



        if len(subdirs) != 0:
            continue

        for fn in files:
            if fn in filelist:
                disk_dict[fn].append(out_subpath)



    for filename in filelist:
        disk_locs = disk_dict[filename]
        sam_locs = samweb.locateFile(filenameorid=filename)
        if len(sam_locs) == 0 and not upload:
            print('No location: %s' % filename)





        locs_to_add = []
        for disk_loc in disk_locs:
            should_add = True
            for sam_loc in sam_locs:
                if sam_loc['location_type'] == 'disk':
                    if disk_loc == sam_loc['location'].split(':')[-1]:
                        should_add = False
                        break
            if should_add:
                locs_to_add.append(disk_loc)







        locs_to_remove = []
        for sam_loc in sam_locs:
            if sam_loc['location_type'] == 'disk':



                if remove:
                    locs_to_remove.append(sam_loc['location'])



                else:



                    local_path = os.path.join(sam_loc['location'].split(':')[-1], filename)
                    if not larbatch_posix.exists(local_path):
                        locs_to_remove.append(sam_loc['location'])






        locs_to_upload = {}
        should_upload = False
        if upload and len(disk_locs) > 0:
            should_upload = True
            for sam_loc in sam_locs:
                if sam_loc['location_type'] == 'tape':
                    should_upload = False
                    break
            if should_upload:
                dropbox = project_utilities.get_dropbox(filename)
                if not larbatch_posix.exists(dropbox):
                    print('Making dropbox directory %s.' % dropbox)
                    larbatch_posix.makedirs(dropbox)
                locs_to_upload[disk_locs[0]] = dropbox



        for loc in locs_to_add:
            node = project_utilities.get_bluearc_server()
            if loc[0:6] == '/pnfs/':
                node = project_utilities.get_dcache_server()
            loc = node + loc.split(':')[-1]
            if add:
                print('Adding location: %s.' % loc)
                project_utilities.test_token()
                samweb.addFileLocation(filenameorid=filename, location=loc)
            elif not upload:
                print('Can add location: %s.' % loc)

        for loc in locs_to_remove:
            if clean or remove:
                print('Removing location: %s.' % loc)
                project_utilities.test_token()
                samweb.removeFileLocation(filenameorid=filename, location=loc)
            elif not upload:
                print('Should remove location: %s.' % loc)

        for loc in list(locs_to_upload.keys()):
            dropbox = locs_to_upload[loc]



            if not larbatch_posix.isdir(dropbox):
                print('Dropbox directory %s does not exist.' % dropbox)
            else:



                dropbox_filename = os.path.join(dropbox, filename)
                if larbatch_posix.exists(dropbox_filename):
                    print('File %s already exists in dropbox %s.' % (filename, dropbox))
                else:



                    loc_filename = os.path.join(loc, filename)
                    print('Copying %s to dropbox directory %s.' % (filename, dropbox))
                    larbatch_posix.copy(loc_filename, dropbox_filename)

    return 0





def docheck_tape(dim):



    result = 0



    import_samweb()



    nbad = 0
    ntot = 0
    filelist = samweb.listFiles(dimensions=dim, stream=True)
    while 1:
        try:
            filename = next(filelist)
        except StopIteration:
            break



        ntot = ntot + 1



        is_on_tape = False
        sam_locs = samweb.locateFile(filenameorid=filename)
        for sam_loc in sam_locs:
            if sam_loc['location_type'] == 'tape':
                is_on_tape = True
                break

        if is_on_tape:
            print('On tape: %s' % filename)
        else:
            result = 1
            nbad = nbad + 1
            print('Not on tape: %s' % filename)

    print('%d files.' % ntot)
    print('%d files need to be store on tape.' % nbad)

    return result





def dojobsub(project, stage, makeup, recur, dryrun, retain):



    jobid = ''



    procmap = ''



    tmpdir = tempfile.mkdtemp()



    tmpworkdir = tempfile.mkdtemp()





    jobsub_workdir_files_args = []



    input_list_name = ''
    if stage.inputlist != '':
        input_list_name = os.path.basename(stage.inputlist)
        work_list_name = os.path.join(tmpworkdir, input_list_name)
        if stage.inputlist != work_list_name:
            input_files = larbatch_posix.readlines(stage.inputlist)
            print('Making input list.')
            work_list = safeopen(work_list_name)
            for input_file in input_files:
                print('Adding input file %s' % input_file)
                work_list.write('%s\n' % input_file.strip())
            work_list.close()
            print('Done making input list.')



    fcls = project.get_fcl(stage.fclname)



    for fcl in fcls:
      if larbatch_posix.exists(fcl):
        workfcl = os.path.join(tmpworkdir, os.path.basename(fcl))
        if os.path.abspath(fcl) != os.path.abspath(workfcl):
          larbatch_posix.copy(fcl, workfcl)







    wrapper_fcl_name = os.path.join(tmpworkdir, 'wrapper.fcl')
    wrapper_fcl = safeopen(wrapper_fcl_name)
    stageNum = 0
    original_project_name = project.name
    original_stage_name = stage.name
    original_project_version = project.version

    for fcl in fcls:
      wrapper_fcl.write('
      wrapper_fcl.write('
      wrapper_fcl.write('\n')



      if stageNum < len(stage.process_name) and stage.process_name[stageNum] != '':
          wrapper_fcl.write('process_name: %s\n' % stage.process_name[stageNum])




      xml_has_metadata = project.file_type != '' or \
                       project.run_type != ''
      if xml_has_metadata:



        if project.release_tag != '':
            wrapper_fcl.write('services.FileCatalogMetadata.applicationVersion: "%s"\n' % \
                                  project.release_tag)
        else:
            wrapper_fcl.write('services.FileCatalogMetadata.applicationVersion: "test"\n')
        if project.file_type:
            wrapper_fcl.write('services.FileCatalogMetadata.fileType: "%s"\n' % \
                              project.file_type)
        if project.run_type:
            wrapper_fcl.write('services.FileCatalogMetadata.runType: "%s"\n' % \
                              project.run_type)




        if stageNum < len(stage.project_name) and stage.project_name[stageNum] != '':
            project.name = stage.project_name[stageNum]
        if stageNum < len(stage.stage_name) and stage.stage_name[stageNum] != '':
            stage.name = stage.stage_name[stageNum]
        if stageNum < len(stage.project_version) and stage.project_version[stageNum] != '':
            project.version = stage.project_version[stageNum]
        sam_metadata = project_utilities.get_sam_metadata(project, stage)
        if sam_metadata:
            wrapper_fcl.write(sam_metadata)
        project.name = original_project_name
        stage.name = original_stage_name
        project.version = original_project_version




      if (not stage.pubs_input and stage.pubs_output) or stage.output_run:
        wrapper_fcl.write('source.firstRun: %d\n' % stage.output_run)





      if stage.maxfluxfilemb != 0 and stageNum == 0:
         wrapper_fcl.write('physics.producers.generator.FluxCopyMethod: "IFDH"\n')
         wrapper_fcl.write('physics.producers.generator.MaxFluxFileMB: %d\n' % stage.maxfluxfilemb)
      wrapper_fcl.write('
      stageNum = 1 + stageNum

    wrapper_fcl.close()






    abssetupscript = project_utilities.get_setup_script_path()
    setupscript = ''
    if not abssetupscript.startswith('/cvmfs/'):
        setupscript = os.path.join(stage.workdir,'setup_experiment.sh')
        larbatch_posix.copy(abssetupscript, setupscript)
        jobsub_workdir_files_args.extend(['-f', setupscript])
        abssetupscript = ''



    if stage.batchname != '':
        workname = stage.batchname
    else:
        workname = '%s-%s-%s' % (stage.name, project.name, project.release_tag)
    workname = workname + os.path.splitext(stage.script)[1]

    workscript = os.path.join(tmpdir, workname)
    if stage.script != workscript:
        larbatch_posix.copy(stage.script, workscript)



    workstartscript = ''
    workstartname = ''
    if stage.start_script != '':
        workstartname = 'start-%s' % workname

        workstartscript = os.path.join(tmpdir, workstartname)
        if stage.start_script != workstartscript:
            larbatch_posix.copy(stage.start_script, workstartscript)



    workstopscript = ''
    workstopname = ''
    if stage.stop_script != '':
        workstopname = 'stop-%s' % workname

        workstopscript = os.path.join(tmpdir, workstopname)
        if stage.stop_script != workstopscript:
            larbatch_posix.copy(stage.stop_script, workstopscript)



    for init_script in stage.init_script:
        if len(init_script) > 0:
            if larbatch_posix.exists(init_script[0]):
                work_init_script = os.path.join(tmpworkdir, os.path.basename(init_script[0]))
                if init_script[0] != work_init_script:
                    larbatch_posix.copy(init_script[0], work_init_script)



    n = len(stage.init_script)
    if n == 0:
        stage.init_script = ''
    else:




        work_init_wrapper = os.path.join(tmpworkdir, 'init_wrapper.sh')
        f = open(work_init_wrapper, 'w')
        f.write('
        for init_script in stage.init_script:
            f.write('echo\n')
            f.write('echo "Executing %s"\n' % os.path.basename(init_script[0]))
            if larbatch_posix.exists(init_script[0]):
                f.write('./%s' % os.path.basename(init_script[0]))
            else:
                f.write(os.path.basename(init_script[0]))
            if len(init_script) > 1:
                f.write(' %s' % ' '.join(init_script[1:]))
            f.write('\n')
            f.write('status=$?\n')
            f.write('echo "%s finished with status $status"\n' % os.path.basename(init_script[0]))
            f.write('if [ $status -ne 0 ]; then\n')
            f.write('  exit $status\n')
            f.write('fi\n')
        f.write('echo\n')
        f.write('echo "Done executing initialization scripts."\n')
        f.close()
        stage.init_script = work_init_wrapper



    for init_source in stage.init_source:
        if init_source != '':
            if larbatch_posix.exists(init_source):
                work_init_source = os.path.join(tmpworkdir, os.path.basename(init_source))
                if init_source != work_init_source:
                    larbatch_posix.copy(init_source, work_init_source)



    n = len(stage.init_source)
    if n == 0:
        stage.init_source = ''
    else:




        work_init_source_wrapper = os.path.join(tmpworkdir, 'init_source_wrapper.sh')
        f = open(work_init_source_wrapper, 'w')
        for init_source in stage.init_source:
            f.write('echo\n')
            f.write('echo "Sourcing %s"\n' % os.path.basename(init_source))
            if larbatch_posix.exists(init_source):
                f.write('source %s\n' % os.path.basename(init_source))
            else:
                f.write('source `which %s`\n' % os.path.basename(init_source))
        f.write('echo\n')
        f.write('echo "Done sourcing initialization scripts."\n')
        f.close()
        stage.init_source = work_init_source_wrapper



    for end_script in stage.end_script:
        if len(end_script) > 0:
            if larbatch_posix.exists(end_script[0]):
                work_end_script = os.path.join(tmpworkdir, os.path.basename(end_script[0]))
                if end_script[0] != work_end_script:
                    larbatch_posix.copy(end_script[0], work_end_script)



    n = len(stage.end_script)
    if n == 0:
        stage.end_script = ''
    else:



        work_end_wrapper = os.path.join(tmpworkdir, 'end_wrapper.sh')
        f = open(work_end_wrapper, 'w')
        f.write('
        for end_script in stage.end_script:
            f.write('echo\n')
            f.write('echo "Executing %s"\n' % os.path.basename(end_script[0]))
            if larbatch_posix.exists(end_script[0]):
                f.write('./%s' % os.path.basename(end_script[0]))
            else:
                f.write(os.path.basename(end_script[0]))
            if len(end_script) > 1:
                f.write(' %s' % ' '.join(end_script[1:]))
            f.write('\n')
            f.write('status=$?\n')
            f.write('echo "%s finished with status $status"\n' % os.path.basename(end_script[0]))
            f.write('if [ $status -ne 0 ]; then\n')
            f.write('  exit $status\n')
            f.write('fi\n')
        f.write('echo\n')
        f.write('echo "Done executing finalization scripts."\n')
        f.close()
        stage.end_script = work_end_wrapper



    for istage in stage.mid_source:
        for mid_source in stage.mid_source[istage]:
            if mid_source != '':
                if larbatch_posix.exists(mid_source):
                    work_mid_source = os.path.join(tmpworkdir, os.path.basename(mid_source))
                    if mid_source != work_mid_source:
                        larbatch_posix.copy(mid_source, work_mid_source)





    if len(stage.mid_source) > 0:
        work_mid_source_wrapper = os.path.join(tmpworkdir, 'mid_source_wrapper.sh')
        f = open(work_mid_source_wrapper, 'w')
        for istage in stage.mid_source:
            for mid_source in stage.mid_source[istage]:
                f.write('if [ $stage -eq %d ]; then\n' % istage)
                f.write('  echo\n')
                f.write('  echo "Sourcing %s"\n' % os.path.basename(mid_source))
                if larbatch_posix.exists(mid_source):
                    f.write('  source %s\n' % os.path.basename(mid_source))
                else:
                    f.write('  source `which %s`\n' % os.path.basename(mid_source))
                f.write('fi\n')
        f.write('echo\n')
        f.write('echo "Done sourcing midstage source initialization scripts for stage $stage."\n')
        f.close()
        stage.mid_source = work_mid_source_wrapper
    else:
        stage.mid_source = ''



    for istage in stage.mid_script:
        for mid_script in stage.mid_script[istage]:
            if len(mid_script) > 0:
                if larbatch_posix.exists(mid_script[0]):
                    work_mid_script = os.path.join(tmpworkdir, os.path.basename(mid_script[0]))
                    if mid_script[0] != work_mid_script:
                        larbatch_posix.copy(mid_script[0], work_mid_script)




    if len(stage.mid_script) > 0:
        work_mid_wrapper = os.path.join(tmpworkdir, 'mid_wrapper.sh')
        f = open(work_mid_wrapper, 'w')
        f.write('
        f.write('stage=$1\n')
        for istage in stage.mid_script:
            for mid_script in stage.mid_script[istage]:
                f.write('if [ $stage -eq %d ]; then\n' % istage)
                f.write('  echo\n')
                f.write('  echo "Executing %s"\n' % os.path.basename(mid_script[0]))
                if larbatch_posix.exists(mid_script[0]):
                    f.write('  ./%s' % os.path.basename(mid_script[0]))
                else:
                    f.write(os.path.basename(mid_script[0]))
                if len(mid_script) > 1:
                    f.write(' %s' % ' '.join(mid_script[1:]))
                f.write('\n')
                f.write('  status=$?\n')
                f.write('  echo "%s finished with status $status"\n' % os.path.basename(mid_script[0]))
                f.write('  if [ $status -ne 0 ]; then\n')
                f.write('    exit $status\n')
                f.write('  fi\n')
                f.write('fi\n')
        f.write('echo\n')
        f.write('echo "Done executing midstage finalization scripts for stage $stage."\n')
        f.close()
        stage.mid_script = work_mid_wrapper
    else:
        stage.mid_script = ''



    helpers = ('root_metadata.py',
               'artroot_filter.py',
               'merge_json.py',
               'subruns.py',
               'validate_in_job.py',
               'mkdir.py',
               'emptydir.py',
               'file_to_url.sh')

    for helper in helpers:



        jobinfo = subprocess.Popen(['which', helper],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        jobout, joberr = jobinfo.communicate()
        jobout = convert_str(jobout)
        joberr = convert_str(joberr)
        rc = jobinfo.poll()
        if rc == 0:
            helper_path = jobout.splitlines()[0].strip()
            work_helper = os.path.join(tmpworkdir, helper)
            if helper_path != work_helper:
                larbatch_posix.copy(helper_path, work_helper)
        else:
            print('Helper script %s not found.' % helper)




    helper_modules = ('larbatch_posix',
                      'project_utilities',
                      'larbatch_utilities',
                      'experiment_utilities',
                      'extractor_dict')

    for helper_module in helper_modules:



        jobinfo = subprocess.Popen(['python'],
                                   stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        cmd = 'import %s\nprint(%s.__file__)\n' % (helper_module, helper_module)
        jobinfo.stdin.write(convert_bytes(cmd))
        jobout, joberr = jobinfo.communicate()
        jobout = convert_str(jobout)
        joberr = convert_str(joberr)
        rc = jobinfo.poll()
        if rc == 0:
            helper_path = jobout.splitlines()[-1].strip()

            work_helper = os.path.join(tmpworkdir, os.path.basename(helper_path))
            if helper_path != work_helper:
                larbatch_posix.copy(helper_path, work_helper)
        else:
            print('Helper python module %s not found.' % helper_module)




    if makeup:

        checked_file = os.path.join(stage.bookdir, 'checked')
        if not larbatch_posix.exists(checked_file):
            raise RuntimeError('Wait for any running jobs to finish and run project.py --check')
        makeup_count = 0



        bad_filename = os.path.join(stage.bookdir, 'bad.list')
        if larbatch_posix.exists(bad_filename):
            lines = larbatch_posix.readlines(bad_filename)
            for line in lines:
                bad_subdir = line.strip()
                if bad_subdir != '':
                    bad_path = os.path.join(stage.outdir, bad_subdir)
                    if larbatch_posix.exists(bad_path):
                        print('Deleting %s' % bad_path)
                        larbatch_posix.rmtree(bad_path)
                    bad_path = os.path.join(stage.logdir, bad_subdir)
                    if larbatch_posix.exists(bad_path):
                        print('Deleting %s' % bad_path)
                        larbatch_posix.rmtree(bad_path)
                    bad_path = os.path.join(stage.bookdir, bad_subdir)
                    if larbatch_posix.exists(bad_path):
                        print('Deleting %s' % bad_path)
                        larbatch_posix.rmtree(bad_path)





        missing_files = []
        if stage.inputdef == '':
            missing_filename = os.path.join(stage.bookdir, 'missing_files.list')
            if larbatch_posix.exists(missing_filename):
                lines = larbatch_posix.readlines(missing_filename)
                for line in lines:
                    words = line.split()
                    if len(words) > 0:
                        missing_files.append(words[0])
            makeup_count = len(missing_files)
            print('Makeup list contains %d files.' % makeup_count)

        if input_list_name != '':
            work_list_name = os.path.join(tmpworkdir, input_list_name)
            if larbatch_posix.exists(work_list_name):
                larbatch_posix.remove(work_list_name)
            work_list = safeopen(work_list_name)
            for missing_file in missing_files:
                work_list.write('%s\n' % missing_file)
            work_list.close()





        if stage.inputdef == '' and stage.inputfile == '' and stage.inputlist == '':
            procs = set(range(stage.num_jobs))




            output_files = os.path.join(stage.bookdir, 'files.list')
            if larbatch_posix.exists(output_files):
                lines = larbatch_posix.readlines(output_files)
                for line in lines:
                    dir = os.path.basename(os.path.dirname(line))
                    dir_parts = dir.split('_')
                    if len(dir_parts) > 1:
                        proc = int(dir_parts[1])
                        if proc in procs:
                            procs.remove(proc)
                if len(procs) != makeup_count:
                    raise RuntimeError('Makeup process list has different length than makeup count.')



                if len(procs) > 0:
                    procmap = 'procmap.txt'
                    procmap_path = os.path.join(tmpworkdir, procmap)
                    procmap_file = safeopen(procmap_path)
                    for proc in procs:
                        procmap_file.write('%d\n' % proc)
                    procmap_file.close()



        import_samweb()



        cpids = []
        cpids_filename = os.path.join(stage.bookdir, 'cpids.list')
        if larbatch_posix.exists(cpids_filename):
            cpids_files = larbatch_posix.readlines(cpids_filename)
            for line in cpids_files:
                cpids.append(line.strip())



        makeup_defname = ''
        if len(cpids) > 0:
            project_utilities.test_token()
            makeup_defname = samweb.makeProjectName(stage.inputdef) + '_makeup'



            cpids_list = ''
            sep = ''
            for cpid in cpids:
                cpids_list = cpids_list + '%s%s' % (sep, cpid)
                sep = ','



            dim = '(defname: %s) minus (consumer_process_id %s and consumed_status consumed)' % (stage.inputdef, cpids_list)



            print('Creating makeup sam dataset definition %s' % makeup_defname)
            project_utilities.test_token()
            samweb.createDefinition(defname=makeup_defname, dims=dim)
            makeup_count = samweb.countFiles(defname=makeup_defname)
            print('Makeup dataset contains %d files.' % makeup_count)



    tarname = 'work%s.tar' % uuid.uuid4()
    tmptar = '%s/%s' % (tmpworkdir, tarname)
    jobinfo = subprocess.Popen(['tar','-cf', tmptar, '-C', tmpworkdir,
                                '--mtime=2018-01-01',
                                '--exclude=%s' % tarname, '.'],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    jobout, joberr = jobinfo.communicate()
    rc = jobinfo.poll()
    if rc != 0:
        raise RuntimeError('Failed to create work tarball in %s' % tmpworkdir)




















    jobsub_workdir_files_args.extend(['--tar_file_name', 'dropbox://%s' % tmptar])






    inputdef = stage.inputdef
    if makeup and makeup_defname != '':
        inputdef = makeup_defname



    prjname = ''
    if inputdef != '':
        import_samweb()
        project_utilities.test_token()
        prjname = samweb.makeProjectName(inputdef)



    mixprjname = ''
    if stage.mixinputdef != '':
        import_samweb()
        project_utilities.test_token()
        mixprjname = 'mix_%s' % samweb.makeProjectName(stage.mixinputdef)



    prj_started = False
    if prjname != '' and stage.prestart != 0:
        ok = project_utilities.start_project(inputdef, prjname,
                                             stage.num_jobs * stage.max_files_per_job,
                                             stage.recur, stage.filelistdef)
        if ok != 0:
            print('Failed to start project.')
            sys.exit(1)
        prj_started = True



    if mixprjname != '' and prj_started:
        ok = project_utilities.start_project(stage.mixinputdef, mixprjname, 0, 0, stage.filelistdef)
        if ok != 0:
            print('Failed to start mix project.')
            sys.exit(1)



    role = project_utilities.get_role()
    if project.role != '':
        role = project.role



    command = ['jobsub_submit']
    command_njobs = 1



    command.append('--group=%s' % project_utilities.get_experiment())
    command.append('--role=%s' % role)
    command.extend(jobsub_workdir_files_args)
    if project.server != '-' and project.server != '':
        command.append('--jobsub-server=%s' % project.server)
    if stage.resource != '':
        command.append('--resource-provides=usage_model=%s' % stage.resource)
    elif project.resource != '':
        command.append('--resource-provides=usage_model=%s' % project.resource)
    if stage.lines != '':
        command.append('--lines=%s' % stage.lines)
    elif project.lines != '':
        command.append('--lines=%s' % project.lines)
    if stage.site != '':
        command.append('--site=%s' % stage.site)
    if stage.blacklist != '':
        command.append('--blacklist=%s' % stage.blacklist)
    if stage.cpu != 0:
        command.append('--cpu=%d' % stage.cpu)
    if stage.disk != '':
        command.append('--disk=%s' % stage.disk)
    if stage.memory != 0:
        command.append('--memory=%d' % stage.memory)
    if project.os != '':
        if stage.singularity == 0:
            command.append('--OS=%s' % project.os)
        else:
            p = project_utilities.get_singularity(project.os)
            if p != '':
                if (stage.num_jobs > 1 or project.force_dag) and \
                   (inputdef != '' or stage.mixinputdef != '') :
                    command.append(r"""--lines='+SingularityImage=\"%s\"'""" % p)
                else:
                    command.append(r"""--lines='+SingularityImage="%s"'""" % p)
            else:
                raise RuntimeError('No singularity image found for %s' % project.os)
    if not stage.pubs_output:
        if not makeup:
            command_njobs = stage.num_jobs
            command.extend(['-N', '%d' % command_njobs])
        else:
            command_njobs = min(makeup_count, stage.num_jobs)
            command.extend(['-N', '%d' % command_njobs])
    else:
        if stage.inputdef != '':
            command_njobs = stage.num_jobs
        else:
            command_njobs = stage.num_jobs
            command.extend(['-N', '%d' % command_njobs])
    if stage.jobsub != '':
        for word in stage.jobsub.split():
            command.append(word)
    opt = project_utilities.default_jobsub_submit_options()
    if opt != '':
        for word in opt.split():
            command.append(word)
    if stage.cvmfs != 0:
        command.append('--append_condor_requirements=\'(TARGET.HAS_CVMFS_%s_opensciencegrid_org==true)\'' % project_utilities.get_experiment())
    if stage.stash != 0:
        command.append('--append_condor_requirements=\'(TARGET.HAS_CVMFS_%s_osgstorage_org==true)\'' % project_utilities.get_experiment())
    if stage.singularity != 0:
        command.append('--append_condor_requirements=\'(TARGET.HAS_SINGULARITY=?=true)\'')



    workurl = "file://%s" % workscript
    command.append(workurl)




    if stage.max_files_per_job != 0:
        command_max_files_per_job = stage.max_files_per_job
        command.extend(['--nfile', '%d' % command_max_files_per_job])




    command.extend([' --group', project_utilities.get_experiment()])
    command.extend([' -g'])
    command.extend([' -c', 'wrapper.fcl'])
    command.extend([' --ups', ','.join(project.ups)])
    if project.release_tag != '':
        command.extend([' -r', project.release_tag])
    command.extend([' -b', project.release_qual])
    if project.local_release_dir != '':
        command.extend([' --localdir', project.local_release_dir])
    if project.local_release_tar != '':
        command.extend([' --localtar', project.local_release_tar])
    command.extend([' --workdir', stage.workdir])
    command.extend([' --outdir', stage.outdir])
    command.extend([' --logdir', stage.logdir])
    if stage.dirsize > 0:
        command.extend([' --dirsize', '%d' % stage.dirsize])
    if stage.dirlevels > 0:
        command.extend([' --dirlevels', '%d' % stage.dirlevels])
    if stage.exe:
        if type(stage.exe) == type([]):
            command.extend([' --exe', ':'.join(stage.exe)])
        else:
            command.extend([' --exe', stage.exe])
    if stage.schema != '':
        command.extend([' --sam_schema', stage.schema])
    if project.os != '':
        command.extend([' --os', project.os])



    if not stage.pubs_input and stage.pubs_output and stage.output_subruns[0] > 0:
        command.extend(['--process', '%d' % (stage.output_subruns[0]-1)])



    if stage.dynamic:
        command.append('--single')

    if stage.inputfile != '':
        command.extend([' -s', stage.inputfile])
    elif input_list_name != '':
        command.extend([' -S', input_list_name])
    elif inputdef != '':
        command.extend([' --sam_defname', inputdef,
                        ' --sam_project', prjname])
    if recur:
        command.extend([' --recur'])
    if stage.mixinputdef != '':
        command.extend([' --mix_defname', stage.mixinputdef,
                        ' --mix_project', mixprjname])
    if stage.inputmode != '':
        command.extend([' --inputmode', stage.inputmode])
    command.extend([' -n', '%d' % stage.num_events])
    if stage.inputdef == '':
        command.extend([' --njobs', '%d' % stage.num_jobs ])
    for ftype in stage.datafiletypes:
        command.extend(['--data_file_type', ftype])
    if procmap != '':
        command.extend([' --procmap', procmap])
    if stage.output:
        if type(stage.output) == type([]):
            command.extend([' --output', ':'.join(stage.output)])
        else:
            command.extend([' --output', stage.output])
    if stage.TFileName != '':
        command.extend([' --TFileName', stage.TFileName])
    if stage.init_script != '':
        command.extend([' --init-script', os.path.basename(stage.init_script)])
    if stage.init_source != '':
        command.extend([' --init-source', os.path.basename(stage.init_source)])
    if stage.end_script != '':
        command.extend([' --end-script', os.path.basename(stage.end_script)])
    if stage.mid_source != '':
        command.extend([' --mid-source', os.path.basename(stage.mid_source)])
    if stage.mid_script != '':
        command.extend([' --mid-script', os.path.basename(stage.mid_script)])
    if abssetupscript != '':
        command.extend([' --init', abssetupscript])



    if stage.validate_on_worker == 1:
      print('Validation will be done on the worker node %d' % stage.validate_on_worker)
      command.extend([' --validate'])
      command.extend([' --declare'])

      if type(stage.fclname) == type([]) and len(stage.fclname) > 1:
        command.extend([' --maintain_parentage'])

    if stage.copy_to_fts == 1:
      command.extend([' --copy'])

    if stage.nthreads > 1:
      command.extend([' --nthreads', '%d' % stage.nthreads])
    if stage.nschedules > 1:
      command.extend([' --nschedules', '%d' % stage.nschedules])



    if (prjname != '' or mixprjname != '') and command_njobs == 1 and not project.force_dag and not prj_started:
        command.extend([' --sam_start',
                        ' --sam_station', project_utilities.get_experiment(),
                        ' --sam_group', project_utilities.get_experiment()])



    if len(stage.args) > 0:
        command.append(' --args')
        command.extend(stage.args)




    start_commands = []
    stop_commands = []
    dag_prjs = []
    if command_njobs > 1 or project.force_dag:
        if inputdef != '':
            dag_prjs.append([inputdef, prjname])
        if stage.mixinputdef != '':
            dag_prjs.append([stage.mixinputdef, mixprjname])

    for dag_prj in dag_prjs:




        if workstartname == '' or workstopname == '':
            raise RuntimeError('Sam start or stop project script not found.')



        start_command = ['jobsub_submit']



        start_command.append('--group=%s' % project_utilities.get_experiment())
        if setupscript != '':
            start_command.append('-f %s' % setupscript)

        if stage.resource != '':
            start_command.append('--resource-provides=usage_model=%s' % stage.resource)
        elif project.resource != '':
            start_command.append('--resource-provides=usage_model=%s' % project.resource)
        if stage.lines != '':
            start_command.append('--lines=%s' % stage.lines)
        elif project.lines != '':
            start_command.append('--lines=%s' % project.lines)
        if stage.site != '':
            start_command.append('--site=%s' % stage.site)
        if stage.blacklist != '':
            start_command.append('--blacklist=%s' % stage.blacklist)
        if project.os != '':
            if stage.singularity == 0:
                start_command.append('--OS=%s' % project.os)
            else:
                p = project_utilities.get_singularity(project.os)
                if p != '':
                    start_command.append('--lines=\'+SingularityImage=\\"%s\\"\'' % p)
                else:
                    raise RuntimeError('No singularity image found for %s' % project.os)
        if stage.jobsub_start != '':
            for word in stage.jobsub_start.split():
                start_command.append(word)
        opt = project_utilities.default_jobsub_submit_options()
        if opt != '':
            for word in opt.split():
                start_command.append(word)
        if stage.cvmfs != 0:
            start_command.append('--append_condor_requirements=\'(TARGET.HAS_CVMFS_%s_opensciencegrid_org==true)\'' % project_utilities.get_experiment())
        if stage.stash != 0:
            start_command.append('--append_condor_requirements=\'(TARGET.HAS_CVMFS_%s_osgstorage_org==true)\'' % project_utilities.get_experiment())
        if stage.singularity != 0:
            start_command.append('--append_condor_requirements=\'(TARGET.HAS_SINGULARITY=?=true)\'')



        workstarturl = "file://%s" % workstartscript
        start_command.append(workstarturl)



        start_command.extend([' --sam_station', project_utilities.get_experiment(),
                              ' --sam_group', project_utilities.get_experiment(),
                              ' --sam_defname', dag_prj[0],
                              ' --sam_project', dag_prj[1],
                              ' -g'])
        if recur:
            start_command.extend([' --recur'])

        if abssetupscript != '':
            start_command.extend([' --init', abssetupscript])

        if stage.num_jobs > 0 and stage.max_files_per_job > 0:
            start_command.extend([' --max_files', '%d' % (stage.num_jobs * stage.max_files_per_job)])

        if stage.prestagefraction > 0.:
            start_command.extend([' --prestage_fraction', '%f' % stage.prestagefraction])



        start_command.extend([' --logdir', stage.logdir])



        if not prj_started or stage.prestagefraction > 0.:
            start_commands.append(start_command)



        stop_command = ['jobsub_submit']



        stop_command.append('--group=%s' % project_utilities.get_experiment())
        if setupscript != '':
            stop_command.append('-f %s' % setupscript)

        if stage.resource != '':
            stop_command.append('--resource-provides=usage_model=%s' % stage.resource)
        elif project.resource != '':
            stop_command.append('--resource-provides=usage_model=%s' % project.resource)
        if stage.lines != '':
            stop_command.append('--lines=%s' % stage.lines)
        elif project.lines != '':
            stop_command.append('--lines=%s' % project.lines)
        if stage.site != '':
            stop_command.append('--site=%s' % stage.site)
        if stage.blacklist != '':
            stop_command.append('--blacklist=%s' % stage.blacklist)
        if project.os != '':
            if stage.singularity == 0:
                stop_command.append('--OS=%s' % project.os)
            else:
                p = project_utilities.get_singularity(project.os)
                if p != '':
                    stop_command.append('--lines=\'+SingularityImage=\\"%s\\"\'' % p)
                else:
                    raise RuntimeError('No singularity image found for %s' % project.os)
        if stage.jobsub_start != '':
            for word in stage.jobsub_start.split():
                stop_command.append(word)
        opt = project_utilities.default_jobsub_submit_options()
        if opt != '':
            for word in opt.split():
                stop_command.append(word)
        if stage.cvmfs != 0:
            stop_command.append('--append_condor_requirements=\'(TARGET.HAS_CVMFS_%s_opensciencegrid_org==true)\'' % project_utilities.get_experiment())
        if stage.stash != 0:
            stop_command.append('--append_condor_requirements=\'(TARGET.HAS_CVMFS_%s_osgstorage_org==true)\'' % project_utilities.get_experiment())
        if stage.singularity != 0:
            stop_command.append('--append_condor_requirements=\'(TARGET.HAS_SINGULARITY=?=true)\'')



        workstopurl = "file://%s" % workstopscript
        stop_command.append(workstopurl)



        stop_command.extend([' --sam_station', project_utilities.get_experiment(),
                             ' --sam_project', dag_prj[1],
                             ' -g'])



        stop_command.extend([' --logdir', stage.logdir])

        if abssetupscript != '':
            stop_command.extend([' --init', abssetupscript])



        stop_commands.append(stop_command)

    if len(start_commands) > 0 or len(stop_commands) > 0:



        dagfilepath = os.path.join(tmpdir, 'submit.dag')
        dag = safeopen(dagfilepath)
        dag.write('<serial>\n')



        if len(start_commands) > 0:
            dag.write('\n<parallel>\n\n')
            for start_command in start_commands:
                first = True
                for word in start_command:
                    if not first:
                        dag.write(' ')
                    dag.write(word)
                    first = False
                dag.write('\n\n')
            dag.write('</parallel>\n')



        dag.write('\n<parallel>\n\n')
        for process in range(command_njobs):

            first = True
            skip = False
            for word in command:
                if skip:
                    skip = False
                else:
                    if word == '-N':

                        skip = True
                    else:
                        if not first:
                            dag.write(' ')
                        if word[:7] == '--role=':
                            word = ''
                        if word.startswith('--jobsub-server='):
                            word = ''
                        word = project_utilities.dollar_escape(word)
                        dag.write(word)
                        first = False
            dag.write(' --process %d\n' % process)
            dag.write('\n')
        dag.write('\n</parallel>\n')



        if len(stop_commands) > 0:
            dag.write('\n<parallel>\n\n')
            for stop_command in stop_commands:
                first = True
                for word in stop_command:
                    if not first:
                        dag.write(' ')
                    dag.write(word)
                    first = False
                dag.write('\n\n')
            dag.write('</parallel>\n')



        dag.write('\n</serial>\n')
        dag.close()



        command = ['jobsub_submit_dag']
        command.append('--group=%s' % project_utilities.get_experiment())
        if project.server != '-' and project.server != '':
            command.append('--jobsub-server=%s' % project.server)
        command.append('--role=%s' % role)
        dagfileurl = 'file://'+ dagfilepath
        command.append(dagfileurl)

    checked_file = os.path.join(stage.bookdir, 'checked')



    submit_timeout = 3600000
    if prjname != '':
        submit_timeout += 1.0 * command_njobs
    if stage.jobsub_timeout > submit_timeout:
        submit_timeout = stage.jobsub_timeout



    if not makeup:



        print('Invoke jobsub_submit')
        if dryrun:
            print(' '.join(command))
        else:
            q = queue.Queue()
            jobinfo = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            thread = threading.Thread(target=project_utilities.wait_for_subprocess, args=[jobinfo, q])
            thread.start()
            thread.join(timeout=submit_timeout)
            if thread.is_alive():
                jobinfo.terminate()
                thread.join()
            rc = q.get()
            jobout = convert_str(q.get())
            joberr = convert_str(q.get())
            if larbatch_posix.exists(checked_file):
                larbatch_posix.remove(checked_file)
            if larbatch_posix.isdir(tmpdir):
                if retain:
                    print('Retaining temporary directory: %s' % tmpdir)
                else:
                    larbatch_posix.rmtree(tmpdir)
            if larbatch_posix.isdir(tmpworkdir):
                if retain:
                    print('Retaining temporary directory: %s' % tmpworkdir)
                else:
                    larbatch_posix.rmtree(tmpworkdir)
            if rc != 0:
                raise JobsubError(command, rc, jobout, joberr)
            for line in jobout.split('\n'):
                if "JobsubJobId" in line:
                    jobid = line.strip().split()[-1]
                elif "Use job id" in line:
                    jobid = line.strip().split()[3]
            print('job id = %s' % jobid)
            if not jobid:
                raise JobsubError(command, rc, jobout, joberr)
        print('jobsub_submit finished.')

    else:



        if makeup_count > 0:
            if dryrun:
                print(' '.join(command))
            else:
                q = queue.Queue()
                jobinfo = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                thread = threading.Thread(target=project_utilities.wait_for_subprocess,
                                          args=[jobinfo, q])
                thread.start()
                thread.join(timeout=submit_timeout)
                if thread.is_alive():
                    jobinfo.terminate()
                    thread.join()
                rc = q.get()
                jobout = convert_str(q.get())
                joberr = convert_str(q.get())
                if larbatch_posix.exists(checked_file):
                    larbatch_posix.remove(checked_file)
                if larbatch_posix.isdir(tmpdir):
                    if retain:
                        print('Retaining temporary directory: %s' % tmpdir)
                    else:
                        larbatch_posix.rmtree(tmpdir)
                if larbatch_posix.isdir(tmpworkdir):
                    if retain:
                        print('Retaining temporary directory: %s' % tmpworkdir)
                    else:
                        larbatch_posix.rmtree(tmpworkdir)
                if rc != 0:
                    raise JobsubError(command, rc, jobout, joberr)
                for line in jobout.split('\n'):
                    if "JobsubJobId" in line:
                        jobid = line.strip().split()[-1]
                if not jobid:
                    raise JobsubError(command, rc, jobout, joberr)
        else:
            print('Makeup action aborted because makeup job count is zero.')



    return jobid




def dosubmit(project, stage, makeup=False, recur=False, dryrun=False, retain=False):



    project_utilities.test_token()



    larbatch_utilities.test_jobsub()



    ok = stage.checksubmit()
    if ok != 0:
        print('No jobs submitted.')
        return





    if stage.pubs_output and not stage.dynamic:
        if larbatch_posix.exists(stage.workdir):
            larbatch_posix.rmtree(stage.workdir)
        if larbatch_posix.exists(stage.outdir):
            larbatch_posix.rmtree(stage.outdir)
        if larbatch_posix.exists(stage.logdir):
            larbatch_posix.rmtree(stage.logdir)
        if larbatch_posix.exists(stage.bookdir):
            larbatch_posix.rmtree(stage.bookdir)



    if not makeup:
        stage.makedirs()
    else:
        stage.checkdirs()



    ok = stage.checkinput(checkdef=True)
    if ok != 0:
        print('No jobs submitted.')
        return



    if not makeup and not recur and not stage.dynamic:
        if len(larbatch_posix.listdir(stage.outdir)) != 0:
            raise RuntimeError('Output directory %s is not empty.' % stage.outdir)
        if len(larbatch_posix.listdir(stage.logdir)) != 0:
            raise RuntimeError('Log directory %s is not empty.' % stage.logdir)
        if len(larbatch_posix.listdir(stage.bookdir)) != 0:
            raise RuntimeError('Log directory %s is not empty.' % stage.bookdir)



    jobid = dojobsub(project, stage, makeup, recur, dryrun, retain)



    jobids_filename = os.path.join(stage.bookdir, 'jobids.list')
    jobids = []
    if larbatch_posix.exists(jobids_filename):
        lines = larbatch_posix.readlines(jobids_filename)
        for line in lines:
            id = line.strip()
            if len(id) > 0:
                jobids.append(id)
    if len(jobid) > 0:
        jobids.append(jobid)

    jobid_file = safeopen(jobids_filename)
    for jobid in jobids:
        jobid_file.write('%s\n' % jobid)
    jobid_file.close()



    return jobid







def domerge(stage, mergehist, mergentuple):

    hlist = []
    hnlist = os.path.join(stage.bookdir, 'filesana.list')
    if larbatch_posix.exists(hnlist):
        hlist = larbatch_posix.readlines(hnlist)
    else:
        raise RuntimeError('No filesana.list file found %s, run project.py --checkana' % hnlist)

    histurlsname_temp = 'histurls.list'
    histurls = safeopen(histurlsname_temp)

    for hist in hlist:
        histurls.write('%s\n' % hist)
    histurls.close()

    if len(hlist) > 0:
        name = os.path.join(stage.outdir, 'anahist.root')
        if name[0:6] == '/pnfs/':
            tempdir = '%s/mergentuple_%d_%d' % (project_utilities.get_scratch_dir(),
                                                os.getuid(),
                                                os.getpid())
            if not larbatch_posix.isdir(tempdir):
                larbatch_posix.makedirs(tempdir)
            name_temp = '%s/anahist.root' % tempdir
        else:
            name_temp = name

        if mergehist:
            mergecom = "hadd -T"
        elif mergentuple:
            mergecom = "hadd"
        else:
            mergecom = stage.merge

        print("Merging %d root files using %s." % (len(hlist), mergecom))

        if larbatch_posix.exists(name_temp):
            larbatch_posix.remove(name_temp)
        comlist = mergecom.split()
        comlist.extend(["-f", "-k", name_temp, '@' + histurlsname_temp])
        rc = subprocess.call(comlist, stdout=sys.stdout, stderr=sys.stderr)
        if rc != 0:
            print("%s exit status %d" % (mergecom, rc))
        if name != name_temp:
            if larbatch_posix.exists(name):
                larbatch_posix.remove(name)
            if larbatch_posix.exists(name_temp):


                larbatch_posix.copy(name_temp, name)
                larbatch_posix.rmtree(tempdir)
        larbatch_posix.remove(histurlsname_temp)




def doaudit(stage):

    import_samweb()
    stage_has_input = stage.inputfile != '' or stage.inputlist != '' or stage.inputdef != ''
    if not stage_has_input:
        raise RuntimeError('No auditing for generator stage.')



    outputlist = []
    outparentlist = []
    if stage.defname != '':
        query = 'isparentof: (defname: %s) and availability: anylocation' %(stage.defname)
        try:
            outparentlist = samweb.listFiles(dimensions=query)
            outputlist = samweb.listFiles(defname=stage.defname)
        except:
            raise RuntimeError('Error accessing sam information for definition %s.\nDoes definition exist?' % stage.defname)
    else:
        raise RuntimeError('Output definition not found.')




    inputlist = []
    if stage.inputdef != '':
        import_samweb()
        inputlist=samweb.listFiles(defname=stage.inputdef)
    elif stage.inputlist != '':
        ilist = []
        if larbatch_posix.exists(stage.inputlist):
            ilist = larbatch_posix.readlines(stage.inputlist)
            inputlist = []
            for i in ilist:
                inputlist.append(os.path.basename(i.strip()))
    else:
        raise RuntimeError('Input definition and/or input list does not exist.')

    difflist = set(inputlist)^set(outparentlist)
    mc = 0;
    me = 0;
    for item in difflist:
        if item in inputlist:
            mc = mc+1
            if mc==1:
                missingfilelistname = os.path.join(stage.bookdir, 'missingfiles.list')
                missingfilelist = safeopen(missingfilelistname)
                if mc>=1:
                    missingfilelist.write("%s\n" %item)
        elif item in outparentlist:
            me = me+1
            childcmd = 'samweb list-files "ischildof: (file_name=%s) and availability: physical"' %(item)
            children = convert_str(subprocess.check_output(childcmd, shell=True)).splitlines()
            rmfile = list(set(children) & set(outputlist))[0]
            if me ==1:
                flist = []
                fnlist = os.path.join(stage.bookdir, 'files.list')
                if larbatch_posix.exists(fnlist):
                    flist = larbatch_posix.readlines(fnlist)
                    slist = []
                    for line in flist:
                        slist.append(line.split()[0])
                else:
                    raise RuntimeError('No files.list file found %s, run project.py --check' % fnlist)



            sdict = {'content_status':'bad'}
            project_utilities.test_token()
            samweb.modifyFileMetadata(rmfile, sdict)
            print('\nDeclaring the status of the following file as bad:', rmfile)



            fn = []
            fn = [x for x in slist if os.path.basename(x.strip()) != rmfile]
            thefile = safeopen(fnlist)
            for item in fn:
                thefile.write("%s\n" % item)

    if mc==0 and me==0:
        print("Everything in order.")
        return 0
    else:
        print('Missing parent file(s) = ', mc)
        print('Extra parent file(s) = ',me)

    if mc != 0:
        missingfilelist.close()
        print("Creating missingfiles.list in the output directory....done!")
    if me != 0:
        thefile.close()

        print("For extra parent files, files.list redefined and content status declared as bad in SAM...done!")




def help():

    filename = sys.argv[0]
    file = open(filename, 'r')

    doprint=0

    for line in file.readlines():
        if line[2:12] == 'project.py':
            doprint = 1
        elif line[0:6] == '
            doprint = 0
        if doprint:
            if len(line) > 2:
                print(line[2:], end=' ')
            else:
                print()













def normxmlpath(xmlfile):



    normxmlfile = xmlfile



    if xmlfile.find(':') < 0 and \
       not xmlfile.startswith('/') and \
       not xmlfile.startswith('./') and \
       not xmlfile.startswith('../') and \
       xmlfile != '-':




        dirs = [os.getcwd()]



        if 'XMLPATH' in os.environ:
            dirs.extend(os.environ['XMLPATH'].split(':'))



        for dir in dirs:
            xmlpath = os.path.join(dir, xmlfile)
            if os.path.exists(xmlpath):
                normxmlfile = xmlpath
                break



    return normxmlfile



def xmlhelp():

    filename = sys.argv[0]
    file = open(filename, 'r')

    doprint=0

    for line in file.readlines():
        if line[2:20] == 'XML file structure':
            doprint = 1
        elif line[0:6] == '
            doprint = 0
        if doprint:
            if len(line) > 2:
                print(line[2:], end=' ')
            else:
                print()




def main(argv):



    xmlfile = ''
    projectname = ''
    stagenames = ['']
    lines = ''
    site = ''
    cpu = 0
    disk = ''
    memory = 0
    inputdef = ''
    merge = 0
    submit = 0
    recur = 0
    pubs = 0
    pubs_run = 0
    pubs_subruns = []
    pubs_version = None
    check = 0
    checkana = 0
    shorten = 0
    fetchlog = 0
    mergehist = 0
    mergentuple = 0
    audit = 0
    retain = 0
    stage_status = 0
    makeup = 0
    clean = 0
    clean_one = 0
    dump_project = 0
    dump_stage = 0
    dryrun = 0
    nocheck = 0
    print_outdir = 0
    print_logdir = 0
    print_workdir = 0
    print_bookdir = 0
    fcl = 0
    defname = 0
    do_input_files = 0
    do_check_submit = 0
    do_check_input = 0
    declare = 0
    declare_ana = 0
    define = 0
    define_ana = 0
    undefine = 0
    check_declarations = 0
    check_declarations_ana = 0
    test_declarations = 0
    test_declarations_ana = 0
    check_definition = 0
    check_definition_ana = 0
    test_definition = 0
    test_definition_ana = 0
    add_locations = 0
    add_locations_ana = 0
    check_locations = 0
    check_locations_ana = 0
    upload = 0
    upload_ana = 0
    check_tape = 0
    check_tape_ana = 0
    clean_locations = 0
    clean_locations_ana = 0
    remove_locations = 0
    remove_locations_ana = 0

    args = argv[1:]
    while len(args) > 0:
        if args[0] == '-h' or args[0] == '--help' :
            help()
            return 0
        elif args[0] == '-xh' or args[0] == '--xmlhelp' :
            xmlhelp()
            return 0
        elif args[0] == '--token':
            del args[0]
        elif args[0] == '--xml' and len(args) > 1:
            xmlfile = args[1]
            del args[0:2]
        elif args[0] == '--project' and len(args) > 1:
            projectname = args[1]
            del args[0:2]
        elif args[0] == '--stage' and len(args) > 1:
            stagenames = args[1].split(',')
            del args[0:2]
        elif args[0] == '--tmpdir' and len(args) > 1:
            os.environ['TMPDIR'] = args[1]
            del args[0:2]
        elif args[0] == '--lines' and len(args) > 1:
            lines = args[1]
            del args[0:2]
        elif args[0] == '--site' and len(args) > 1:
            site = args[1]
            del args[0:2]
        elif args[0] == '--cpu' and len(args) > 1:
            cpu = int(args[1])
            del args[0:2]
        elif args[0] == '--disk' and len(args) > 1:
            disk = args[1]
            del args[0:2]
        elif args[0] == '--memory' and len(args) > 1:
            memory = int(args[1])
            del args[0:2]
        elif args[0] == '--inputdef' and len(args) > 1:
            inputdef = args[1]
            del args[0:2]
        elif args[0] == '--submit':
            submit = 1
            del args[0]
        elif args[0] == '--recur':
            recur = 1
            del args[0]
        elif args[0] == '--pubs' and len(args) > 2:
            pubs = 1
            pubs_run = int(args[1])
            pubs_subruns = project_utilities.parseInt(args[2])
            del args[0:3]
            if len(args) > 0 and args[0] != '' and args[0][0] != '-':
                pubs_version = int(args[0])
                del args[0]
        elif args[0] == '--check':
            check = 1
            del args[0]
        elif args[0] == '--checkana':
            checkana = 1
            del args[0]
        elif args[0] == '--shorten':
            shorten = 1
            del args[0]
        elif args[0] == '--fetchlog':
            fetchlog = 1
            del args[0]
        elif args[0] == '--merge':
            merge = 1
            del args[0]
        elif args[0] == '--mergehist':
            mergehist = 1
            del args[0]
        elif args[0] == '--mergentuple':
            mergentuple = 1
            del args[0]
        elif args[0] == '--audit':
            audit = 1
            del args[0]
        elif args[0] == '--retain':
            retain = 1
            del args[0]
        elif args[0] == '--status':
            stage_status = 1
            del args[0]
        elif args[0] == '--makeup':
            makeup = 1
            del args[0]
        elif args[0] == '--clean':
            clean = 1
            del args[0]
        elif args[0] == '--clean_one':
            clean_one = 1
            del args[0]
        elif args[0] == '--dump_project':
            dump_project = 1
            del args[0]
        elif args[0] == '--dump_stage':
            dump_stage = 1
            del args[0]
        elif args[0] == '--dryrun':
            dryrun = 1
            del args[0]
        elif args[0] == '--nocheck':
            nocheck = 1
            del args[0]
        elif args[0] == '--outdir':
            print_outdir = 1
            del args[0]
        elif args[0] == '--logdir':
            print_logdir = 1
            del args[0]
        elif args[0] == '--workdir':
            print_workdir = 1
            del args[0]
        elif args[0] == '--bookdir':
            print_bookdir = 1
            del args[0]
        elif args[0] == '--fcl':
            fcl = 1
            del args[0]
        elif args[0] == '--defname':
            defname = 1
            del args[0]
        elif args[0] == '--input_files':
            do_input_files = 1
            del args[0]
        elif args[0] == '--check_submit':
            do_check_submit = 1
            del args[0]
        elif args[0] == '--check_input':
            do_check_input = 1
            del args[0]
        elif args[0] == '--declare':
            declare = 1
            del args[0]
        elif args[0] == '--declare_ana':
            declare_ana = 1
            del args[0]
        elif args[0] == '--define':
            define = 1
            del args[0]
        elif args[0] == '--define_ana':
            define_ana = 1
            del args[0]
        elif args[0] == '--undefine':
            undefine = 1
            del args[0]
        elif args[0] == '--check_declarations':
            check_declarations = 1
            del args[0]
        elif args[0] == '--check_declarations_ana':
            check_declarations_ana = 1
            del args[0]
        elif args[0] == '--test_declarations':
            test_declarations = 1
            del args[0]
        elif args[0] == '--test_declarations_ana':
            test_declarations_ana = 1
            del args[0]
        elif args[0] == '--check_definition':
            check_definition = 1
            del args[0]
        elif args[0] == '--check_definition_ana':
            check_definition_ana = 1
            del args[0]
        elif args[0] == '--test_definition':
            test_definition = 1
            del args[0]
        elif args[0] == '--test_definition_ana':
            test_definition_ana = 1
            del args[0]
        elif args[0] == '--add_locations':
            add_locations = 1
            del args[0]
        elif args[0] == '--add_locations_ana':
            add_locations_ana = 1
            del args[0]
        elif args[0] == '--check_locations':
            check_locations = 1
            del args[0]
        elif args[0] == '--check_locations_ana':
            check_locations_ana = 1
            del args[0]
        elif args[0] == '--upload':
            upload = 1
            del args[0]
        elif args[0] == '--upload_ana':
            upload_ana = 1
            del args[0]
        elif args[0] == '--check_tape':
            check_tape = 1
            del args[0]
        elif args[0] == '--check_tape_ana':
            check_tape_ana = 1
            del args[0]
        elif args[0] == '--clean_locations':
            clean_locations = 1
            del args[0]
        elif args[0] == '--clean_locations_ana':
            clean_locations_ana = 1
            del args[0]
        elif args[0] == '--remove_locations':
            remove_locations = 1
            del args[0]
        elif args[0] == '--remove_locations_ana':
            remove_locations_ana = 1
            del args[0]
        else:
            print('Unknown option %s' % args[0])
            return 1



    project_utilities.use_token_auth()



    xmlfile = normxmlpath(xmlfile)



    if xmlfile == '':
        print('No xml file specified.  Type "project.py -h" for help.')
        return 1




    num_action = submit + check + checkana + fetchlog + merge + mergehist + mergentuple + audit + stage_status + makeup + define + define_ana + undefine + declare + declare_ana
    if num_action > 1:
        print('More than one action was specified.')
        return 1



    projects = get_projects(xmlfile, check=(not nocheck))



    for stagename in stagenames:
        project = select_project(projects, projectname, stagename)
        if project != None:
            if projectname == '':
                projectname = project.name
        else:
            raise RuntimeError('No project selected.\n')



    if clean:
        for stagename in stagenames:
            docleanx(projects, projectname, stagename, clean_descendants = True)



    if clean_one:
        for stagename in stagenames:
            docleanx(projects, projectname, stagename, clean_descendants = False)



    if stage_status:
        dostatus(projects)
        return 0




    stages = {}
    for stagename in stagenames:
        stage = project.get_stage(stagename)
        stages[stagename] = stage



        if lines != '':
            stage.lines = lines
        if site != '':
            stage.site = site
        if cpu != 0:
            stage.cpu = cpu
        if disk != '':
            stage.disk = disk
        if memory != 0:
            stage.memory = memory
        if inputdef != '':
            stage.inputdef = inputdef
            stage.inputfile = ''
            stage.inputlist = ''
        if recur != 0:
            stage.recur = recur



        if pubs:
            stage.pubsify_input(pubs_run, pubs_subruns, pubs_version)
            stage.pubsify_output(pubs_run, pubs_subruns, pubs_version)



        if stage.recur and stage.inputdef != '' and stage.basedef != '':



            import_samweb()
            def_exists = False
            try:
                desc = samweb.descDefinition(defname=stage.inputdef)
                def_exists = True
            except samweb_cli.exceptions.DefinitionNotFound:
                pass

            if not def_exists:



                project_utilities.test_token()



                dim = ''



                project_wildcard = '%s_%%' % samweb.makeProjectName(stage.inputdef).rsplit('_',1)[0]
                if stage.recurtype == 'snapshot':
                    dim = 'defname: %s minus snapshot_for_project_name %s' % \
                        (stage.basedef, project_wildcard)
                elif stage.recurtype == 'consumed':
                    dim = 'defname: %s minus (project_name %s and consumed_status consumed)' % \
                        (stage.basedef, project_wildcard)

                elif stage.recurtype == 'child':




                    nstream = 1
                    if stage.data_stream != None and len(stage.data_stream) > 0:
                        nstream = len(stage.data_stream)

                    dim = ''
                    for istream in range(nstream):
                        idim = project_utilities.dimensions_datastream(project, stage,
                                                                       ana=False, index=istream)
                        if idim.find('anylocation') > 0:
                            idim = idim.replace('anylocation', 'physical')
                        else:
                            idim += ' with availability physical'

                        if len(dim) > 0:
                            dim += ' or '
                        dim += '(defname: %s minus isparentof:( %s ) )' % (stage.basedef, idim)

                    if stage.activebase != '':
                        activedef = '%s_active' % stage.activebase
                        waitdef = '%s_wait' % stage.activebase
                        dim += ' minus defname: %s' % activedef
                        dim += ' minus defname: %s' % waitdef
                        project_utilities.makeDummyDef(activedef)
                        project_utilities.makeDummyDef(waitdef)

                elif stage.recurtype == 'anachild':




                    nstream = 1
                    if stage.ana_data_stream != None and len(stage.ana_data_stream) > 0:
                        nstream = len(stage.ana_data_stream)

                    dim = ''
                    for istream in range(nstream):
                        idim = project_utilities.dimensions_datastream(project, stage,
                                                                       ana=True, index=istream)
                        if idim.find('anylocation') > 0:
                            idim = idim.replace('anylocation', 'physical')
                        else:
                            idim += ' with availability physical'

                        if len(dim) > 0:
                            dim += ' or '
                        dim += '(defname: %s minus isparentof:( %s ) )' % (stage.basedef, idim)

                    if stage.activebase != '':
                        activedef = '%s_active' % stage.activebase
                        waitdef = '%s_wait' % stage.activebase
                        dim += ' minus defname: %s' % activedef
                        dim += ' minus defname: %s' % waitdef
                        project_utilities.makeDummyDef(activedef)
                        project_utilities.makeDummyDef(waitdef)

                elif stage.recurtype != '' and stage.recurtype != 'none':
                    raise RuntimeError('Unknown recursive type %s.' % stage.recurtype)



                if stage.recurlimit != 0:
                    dim += ' with limit %d' % stage.recurlimit



                print('Creating recursive dataset definition %s' % stage.inputdef)
                project_utilities.test_token()
                samweb.createDefinition(defname=stage.inputdef, dims=dim)



        valok = larbatch_utilities.validate_stage(project, stage)
        if valok:
            print('Validation checks passed.')
        else:
            print('Validation checks failed.')
            print('Qutting now.')
            sys.exit(1)



    if dump_stage:
        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            print(stage)



    if dump_project:
        print(project)



    if print_outdir:
        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            print(stage.outdir)



    if print_logdir:
        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            print(stage.logdir)



    if print_workdir:
        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            print(stage.workdir)



    if print_bookdir:
        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            print(stage.bookdir)



    if defname:
        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            if stage.defname != '':
                print(stage.defname)



    if do_input_files:
        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            input_files = get_input_files(stage)
            for input_file in input_files:
                print(input_file)



    if do_check_submit:
        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            stage.checksubmit()



    if do_check_input:
        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            stage.checkinput(checkdef=True)



    if shorten:
        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            doshorten(stage)



    rc = 0

    if submit or makeup:



        for stagename in stagenames:
            print('Stage %s:' % stagename)

            if project_utilities.check_running(xmlfile, stagename):
                print('Skipping job submission because similar job submission process is running.')
            else:
                stage = stages[stagename]
                dosubmit(project, stage, makeup, stage.recur, dryrun, retain)

    if check or checkana:



        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            docheck(project, stage, checkana or stage.ana, stage.validate_on_worker)

    if fetchlog:



        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            rc += dofetchlog(project, stage)

    if mergehist or mergentuple or merge:




        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            domerge(stage, mergehist, mergentuple)

    if audit:



        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            doaudit(stage)

    if check_definition or define:



        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            if stage.ana:
                if stage.ana_defname == '':
                    print('No sam analysis dataset definition name specified for this stage.')
                    return 1
                dim = project_utilities.dimensions_datastream(project, stage, ana=True)
                docheck_definition(stage.ana_defname, dim, define)
            else:
                if stage.defname == '':
                    print('No sam dataset definition name specified for this stage.')
                    return 1
                dim = project_utilities.dimensions_datastream(project, stage, ana=False)
                docheck_definition(stage.defname, dim, define)

    if check_definition_ana or define_ana:



        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            if stage.ana_defname == '':
                print('No sam analysis dataset definition name specified for this stage.')
                return 1
            dim = project_utilities.dimensions_datastream(project, stage, ana=True)
            docheck_definition(stage.ana_defname, dim, define_ana)

    if test_definition:



        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            if stage.ana:
                if stage.ana_defname == '':
                    print('No sam dataset definition name specified for this stage.')
                    return 1
                rc += dotest_definition(stage.ana_defname)
            else:
                if stage.defname == '':
                    print('No sam dataset definition name specified for this stage.')
                    return 1
                rc += dotest_definition(stage.defname)

    if test_definition_ana:



        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            if stage.ana_defname == '':
                print('No sam dataset definition name specified for this stage.')
                return 1
            rc += dotest_definition(stage.ana_defname)

    if undefine:



        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            if stage.defname == '':
                print('No sam dataset definition name specified for this stage.')
                return 1
            rc += doundefine(stage.defname)

    if check_declarations or declare:



        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            docheck_declarations(stage.bookdir, stage.outdir, declare, ana=stage.ana)

    if check_declarations_ana or declare_ana:



        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            docheck_declarations(stage.bookdir, stage.outdir, declare_ana, ana=True)

    if test_declarations:



        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            dim = project_utilities.dimensions_datastream(project, stage, ana=stage.ana)
            rc += dotest_declarations(dim)

    if test_declarations_ana:



        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            dim = project_utilities.dimensions_datastream(project, stage, ana=True)
            rc += dotest_declarations(dim)

    if check_locations or add_locations or clean_locations or remove_locations or upload:



        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            dim = project_utilities.dimensions_datastream(project, stage, ana=stage.ana)
            docheck_locations(dim, stage.outdir,
                              add_locations, clean_locations, remove_locations,
                              upload)

    if check_locations_ana or add_locations_ana or clean_locations_ana or \
       remove_locations_ana or upload_ana:



        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            dim = project_utilities.dimensions_datastream(project, stage, ana=True)
            docheck_locations(dim, stage.outdir,
                              add_locations_ana, clean_locations_ana, remove_locations_ana,
                              upload_ana)

    if check_tape:



        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            dim = project_utilities.dimensions_datastream(project, stage, ana=stage.ana)
            docheck_tape(dim)

    if check_tape_ana:



        for stagename in stagenames:
            print('Stage %s:' % stagename)
            stage = stages[stagename]
            dim = project_utilities.dimensions_datastream(project, stage, ana=True)
            docheck_tape(dim)



    return rc



def safeopen(destination):
    if larbatch_posix.exists(destination):
        larbatch_posix.remove(destination)
    file = larbatch_posix.open(destination, 'w')
    return file




def scan_file(fileName):

    returnArray = []
    try:

        fileList = project_utilities.saferead(fileName)

    except:

        return [ -1 ]

    if len(fileList) > 0:
        for line in fileList:
            returnArray.append(line.strip())

    else:


        return [ -1 ]

    return returnArray

if __name__ == '__main__':
    sys.exit(main(sys.argv))

'''inputlist = []
                inp = open(stage.inputlist,"r")
                for line in inp:
                    columns = line.split("/")
                    columns = [col.strip() for col in columns]
                    inputlist.append(columns[8])
                inp.close()'''