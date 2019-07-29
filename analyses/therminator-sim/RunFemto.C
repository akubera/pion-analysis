///
/// \file RunFemto.C
///

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
// #include <TTree.h>

#include "ThermSim.hh"

#include <mutex>
#include <thread>
#include <iostream>
#include <optional>
#include <random>



bool select_particle(TLorentzVector x, TLorentzVector p)
{
  const double
    pt = p.Pt(),
    eta = p.Eta();

  if (!(0.2 < pt && pt < 2.0)) {
    return false;
  }
  if (!(-0.8 < eta && eta < 0.8)) {
    return false;
  }

  // if (0.2 <= x.Perp2() || 0.2 <= x.Z()) {
  //   return false;
  // }

  return true;
}


void
RunFemto(TTree &tree)
{
  std::cout << tree.GetName() << "\n";

  // TTreeReader rg(tree);
  // TTreeReaderValue<float> px(rg, "px");
  TTreeReader tr(&tree);
  TTreeReaderValue<float>
    E(tr, "particle.e"),
    px(tr, "particle.px"),
    py(tr, "particle.py"),
    pz(tr, "particle.pz"),

    t(tr, "particle.t"),
    x(tr, "particle.x"),
    y(tr, "particle.y"),
    z(tr, "particle.z");

  TTreeReaderValue<UInt_t>
    eid(tr, "particle.eventid");

  TTreeReaderValue<Int_t>
    pid(tr, "particle.pid");

  // TTreeReaderValue<Particle>
  //   part(tr, "particle");

  UInt_t prev_eid = -1;
  std::deque<std::unique_ptr<Event>> event_queue;

  std::vector<Particle> particle_vector;
  std::mutex queue_mutex;
  std::condition_variable cv;

  const UInt_t N = 8;

  int elimit = -1;

  std::vector<FemtoRunner> runners;
  runners.reserve(N);
  for (UInt_t i=0; i<N; ++i) {
    runners.emplace_back(event_queue, queue_mutex, cv);
  }

  std::vector<std::thread> threads;
  threads.reserve(N);
  for (auto &runner : runners) {
    threads.emplace_back(runner.spawn());
  }

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> get_rndm_phi(0.0, 2.0 * M_PI);

  double phi = NAN;
  while (tr.Next()) {

    if (prev_eid != *eid) {
        // event_queue.({});
        if (!particle_vector.empty()) {
          std::unique_lock<std::mutex> _lock(queue_mutex);
          auto ev = std::make_unique<Event>(std::move(particle_vector));

          event_queue.emplace_back(std::move(ev));
          cv.notify_one();
        }
        // assert(particle_vector.size() == 0);

        if (elimit-- == 0) {
          break;
        }

        phi = get_rndm_phi(gen);
        prev_eid = *eid;
    }

    TLorentzVector
      vx(*x, *y, *z, *t),
      vp(*px, *py, *pz, *E);

    vx.RotateZ(phi);
    vp.RotateZ(phi);

    if (*pid != 211) {
      continue;
    }

    if (select_particle(vx, vp)) {
      particle_vector.emplace_back(vx, vp);
    }
  }

  for (UInt_t i=0; i<N; ++i) {
    event_queue.emplace_back(nullptr);
  }

  std::cout << "\n";
  while (event_queue.size()) {
    std::cout << "\r" << event_queue.size() - N;
    usleep(300);
  }
  std::cout << "\n";

  for (auto &thread : threads) {
    thread.join();
  }

  for (size_t i=1; i<N; ++i) {
    runners[0].num->Add(runners[i].num.get());
    runners[0].den->Add(runners[i].den.get());
  }

  TFile ofile("ofile.root", "RECREATE");

  runners[0].num->Write();
  runners[0].den->Write();
}
